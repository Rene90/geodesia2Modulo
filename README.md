# 📝 **Programa de Geodesia Geométrica en Rust**

Este programa resuelve problemas fundamentales de geodesia geométrica, incluyendo transformaciones de coordenadas y cálculos de distancias/azimuts usando los métodos de **Puissant** (para distancias cortas) y **Bessel** (para precisión elipsoidal en distancias largas).

## 📋 **Características Principales**
1. **Transformación de Coordenadas:**
   - Geodésicas (latitud, longitud, altura) ↔ Cartesianas (X, Y, Z)
2. **Problemas Geodésicos:**
   - **Directo:** Calcular coordenadas destino desde un punto + azimut/distancia
   - **Inverso:** Calcular distancia y azimut entre dos puntos
3. **Métodos Implementados:**
   - Fórmula de Puissant (para <150 km)
   - Método de Bessel (solución exacta para el elipsoide)
4. **Cálculos Auxiliares:**
   - Radios de curvatura (meridiano `M` y primer vertical `N`)
   - Longitud de arcos (meridiano y paralelo)
   - Aproximación a una esfera

---

## 🛠 **Requisitos**
- **Rust** (instalado via [rustup](https://rustup.rs/))
- Dependencias (en `Cargo.toml`):
  ```toml
  [dependencies]
  serde = { version = "1.0", features = ["derive"] }
  csv = "1.2"

  🚀 Cómo Usar
1. Ejecución Básica
bash
cargo run
El programa pedirá:

Coordenadas de dos puntos en grados decimales (latitud, longitud, altura).

2. Ejemplo de Entrada/Salida
plaintext
Ingrese las coordenadas del primer punto:
Latitud: 19.4326°  
Longitud: -99.1332°  
Altura: 2240 m

Ingrese las coordenadas del segundo punto:
Latitud: 20.6736°  
Longitud: -103.344°  
Altura: 1580 m

Resultados:
- Distancia (Puissant): 450.25 km
- Azimut (Puissant): 45.78°
- Distancia (Bessel): 450.30 km
- Azimut (Bessel): 45.75°

3. Desde Archivo CSV
Modifica la ruta en leer_puntos_archivo() para procesar datos de GPS (formato: Longitude, Latitude, Ellipsoidal height).

🤝 Contribuir
Abre un  issue o envía un pull request!

Código basado en las Notas de Lectura 26 del departamento de Geomatica y Geodesia de la Universidad de New Brunswick (Krakiwsky, E. J., & Thomson, D. B. (1974). *Geodetic Position Computations* (3rd ed.). Department of Surveying Engineering, University of New Brunswick. https://www2.unb.ca/gge/Pubs/LN26.pdf) y en el libro de  Vaníček, P., & Krakiwsky, E. J. (1986). *Geodesy: The Concepts* (2nd ed.). Elsevier Science. https://doi.org/10.1016/B978-0-444-87775-5.X5001-4