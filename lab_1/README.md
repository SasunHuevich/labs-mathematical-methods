# Лабораторная 1 — размещение VM (Julia)

- **[README_DEFENSE.md](README_DEFENSE.md)** — **шпаргалка для защиты** (что говорить, вопросы и ответы).
- **[README_RESULTS.md](README_RESULTS.md)** — результаты запуска (Z, статусы, таблицы).
- **[README_EXPLANATION.md](README_EXPLANATION.md)** — подробно: постановка, код, установка, отладка.

## Быстрый старт

```powershell
cd lab_1
julia install_pkgs.jl
julia run_lab1.jl
```

Блокнот: `Lab1_VM_Placement.ipynb` (ядро **Julia (lab1)** после `IJulia.installkernel(...)`).

## Главный результат

Эталонная модель `data/model.mps` решается HiGHS; **Z ≈ 203.35** за 600 с (допустимое решение, возможен более низкий Z при большем лимите времени).
