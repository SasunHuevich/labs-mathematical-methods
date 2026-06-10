using JuMP
using HiGHS

const DATA_DIR = abspath(joinpath(@__DIR__, "data"))

function parse_general_info(path)
    text = read(path, String)
    params = Dict{String, Any}()
    for m in eachmatch(r"([^=,]+)\s*=\s*([\d.]+)", text)
        params[strip(m.captures[1])] = parse(Float64, m.captures[2])
    end
    params["n"] = Int(params["n"])
    params["m"] = Int(params["m"])
    params["G"] = Int(params["|G|"])
    params["locVMs"] = Int(params["locVMs"])
    params["R"] = Int(params["|R|"])
    params["V"] = Int(params["|V|"])
    params["eps"] = Int(params["eps"])

    plan_vals = Float64[]
    in_plan = false
    for line in split(text, '\n')
        s = strip(line)
        if occursin("Plan_Gg", s)
            in_plan = true
            continue
        end
        in_plan || continue
        isempty(s) && continue
        s == "]" && break
        try
            push!(plan_vals, parse(Float64, replace(s, "[" => "", "]" => "")))
        catch
            break
        end
    end
    params["Plan_G"] = plan_vals[1:Int(params["G"])]
    params["conflict_pairs"] = [(4, 5), (6, 7), (8, 9)]
    return params
end

function read_matrix(path)
    rows = Vector{Vector{Float64}}()
    for line in eachline(path)
        s = strip(line)
        isempty(s) && continue
        push!(rows, parse.(Float64, split(s)))
    end
    return hcat(rows...)'
end

function read_J_for_G(path)
    text = read(path, String)
    i = findfirst("[[", text)
    j = findlast("]]", text)
    inner = text[i[2]+1:j[1]-1]
    servers = Vector{Int}[]
    for chunk in split(inner, "], [")
        chunk = strip(replace(chunk, "[" => "", "]" => ""))
        isempty(chunk) && continue
        push!(servers, parse.(Int, split(chunk, ",")))
    end
    return servers
end

function read_initial_locations(path)
    locs = Tuple{Int, Int}[]
    for line in eachline(path)
        s = strip(line)
        occursin(r"^\(\d+\s*,\s*\d+\)", s) || continue
        m = match(r"\((\d+)\s*,\s*(\d+)\)", s)
        push!(locs, (parse(Int, m.captures[1]), parse(Int, m.captures[2])))
    end
    return locs
end

function load_data(data_dir=DATA_DIR)
    info = parse_general_info(joinpath(data_dir, "general_info.txt"))
    c_j = read_matrix(joinpath(data_dir, "c_y_j.txt"))[:]
    c_ij = read_matrix(joinpath(data_dir, "C_x_ij.txt"))
    d = read_matrix(joinpath(data_dir, "D_ir.txt"))
    q = read_matrix(joinpath(data_dir, "Q_jr.txt"))
    J_for_G = read_J_for_G(joinpath(data_dir, "J_for_G.txt"))
    initial_locs = read_initial_locations(joinpath(data_dir, "x_ij_firstlyLocated.txt"))

    n, m = info["n"], info["m"]
    @assert size(c_ij) == (n, m)
    @assert size(d) == (n, info["R"])
    @assert size(q) == (m, info["R"])
    @assert length(c_j) == m
    @assert length(J_for_G) == info["G"]

    return (
        n = n,
        m = m,
        G = info["G"],
        R = info["R"],
        alpha = info["alpha"],
        beta = info["beta"],
        eps = info["eps"],
        Plan_G = info["Plan_G"],
        conflict_pairs = info["conflict_pairs"],
        c_j = c_j,
        c_ij = c_ij,
        d = d,
        q = q,
        J_for_G = J_for_G,
        initial_locs = initial_locs,
    )
end

function build_vm_model(data; time_limit=3600.0)
    n, m = data.n, data.m
    alpha, beta, eps = data.alpha, data.beta, data.eps

    model = Model(HiGHS.Optimizer)
    set_attribute(model, "time_limit", time_limit)
    set_silent(model)

    allowed = [(i, j) for i in 1:n, j in 1:m if data.c_ij[i, j] > 0]
    allowed_j = [Int[] for _ in 1:n]
    allowed_i = [Int[] for _ in 1:m]
    for (i, j) in allowed
        push!(allowed_j[i], j)
        push!(allowed_i[j], i)
    end

    @variable(model, y[1:m], Bin)
    @variable(model, x[i=1:n, j=1:m], Bin)

    for i in 1:n, j in 1:m
        if data.c_ij[i, j] == 0
            @constraint(model, x[i, j] == 0)
        end
    end

    @objective(
        model,
        Min,
        alpha * sum(data.c_j[j] * y[j] for j in 1:m) +
            beta * sum(data.c_ij[i, j] * x[i, j] for (i, j) in allowed),
    )

    for i in 1:n
        @constraint(model, sum(x[i, j] for j in allowed_j[i]) == 1)
    end

    for j in 1:m, r in 1:data.R
        @constraint(
            model,
            sum(data.d[i, r] * x[i, j] for i in allowed_i[j]) <= data.q[j, r] * y[j],
        )
    end

    for (i1, i2) in data.conflict_pairs, j in 1:m
        if data.c_ij[i1, j] > 0 && data.c_ij[i2, j] > 0
            @constraint(model, x[i1, j] + x[i2, j] <= 1)
        end
    end

    for (g, servers) in enumerate(data.J_for_G)
        @constraint(
            model,
            sum(data.c_j[j] * y[j] for j in servers) +
            sum(data.c_ij[i, j] * x[i, j] for i in 1:n, j in servers) == data.Plan_G[g],
        )
    end

    @constraint(
        model,
        sum(x[i, j] for (i, j) in data.initial_locs if data.c_ij[i, j] > 0) >= eps,
    )

    return model
end

function extract_solution(model, data)
    n, m = data.n, data.m
    placement = Dict{Int, Int}()
    for i in 1:n, j in 1:m
        if data.c_ij[i, j] > 0 && value(model[:x][i, j]) > 0.5
            placement[i] = j
        end
    end
    active_servers = [j for j in 1:m if value(model[:y][j]) > 0.5]
    stayed = count(((i, j),) -> get(placement, i, -1) == j, data.initial_locs)
    return placement, active_servers, stayed
end

function main()
    println("=== Загрузка данных ===")
    data = load_data()
    println("OK: n=$(data.n), m=$(data.m), |G|=$(data.G), locs=$(length(data.initial_locs))")
    allowed = count(data.c_ij .> 0)
    println("Разрешённых пар (i,j): $allowed")

    println("\n=== Построение модели ===")
    t_build = @elapsed begin
        global model = build_vm_model(data; time_limit=600.0)
    end
    println("OK: переменных=$(num_variables(model)), время=$(round(t_build, digits=1)) с")

    println("\n=== Решение HiGHS (лимит 600 с) ===")
    t_solve = @elapsed optimize!(model)
    term = termination_status(model)
    println("Статус: $term, время=$(round(t_solve, digits=1)) с")

    if has_values(model)
        println("Z = ", objective_value(model))
        placement, active_servers, stayed = extract_solution(model, data)
        println("Активных серверов: $(length(active_servers)) / $(data.m)")
        println("VM на исходных местах: $stayed (>= $(data.eps))")
    else
        println("Оптимальное решение за 600 с не получено.")
    end

    println("\n=== Проверка MPS (лимит 600 с) ===")
    try
        t_mps = @elapsed begin
            global mps_model = read_from_file(joinpath(DATA_DIR, "model.mps"))
            set_optimizer(mps_model, HiGHS.Optimizer)
            set_attribute(mps_model, "time_limit", 600.0)
            set_silent(mps_model)
            optimize!(mps_model)
        end
        println("MPS статус: $(termination_status(mps_model)), время=$(round(t_mps, digits=1)) с")
        has_values(mps_model) && println("MPS objective = $(objective_value(mps_model))")
    catch e
        println("MPS пропущен: ", sprint(showerror, e))
    end
end

if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    main()
end
