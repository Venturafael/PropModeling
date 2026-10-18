# PropModeling

Repositório acadêmico desenvolvido a partir do TCC **“Modelagem de Hélices: da comparação entre métodos clássicos ao projeto ótimo”**.

Reúne implementações em **MATLAB e Python** dos métodos BET, BEMT e Vortex Theory para análise de desempenho de hélices, além do método **Minimum Energy Loss (MEL)** de Adkins e Liebeck para projeto ótimo.

Os códigos utilizam polares Clark–Y e dados de referência da APC para calcular e comparar tração, potência e eficiência, incluindo análises de erro e tempo de processamento.

## Estrutura

- `PropModelingMatlab/` — códigos MATLAB.
- `PropModelingPy/` — códigos Python.
- `TCC_Rafael_Ventura.pdf` — fundamentação teórica e resultados.

## Execução

**MATLAB**
```matlab
CompareMethods
```

**Python**
```bash
python3 compare_methods.py
python3 compare_all_apc.py --jobs 4
python3 design_optimal_propeller.py
```
