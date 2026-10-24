# Arquitectura S3 de vectorialcalculus

**Estado:** propuesta para revisión (PR-00; no implementada).  
**Referencia:** vectorialcalculus 1.0.5, rama `master`, 24 de octubre de 2026.  
**Alcance:** diseño de objetos, contratos matemáticos, regiones, evaluación, discretización y migración compatible.  
**No implica:** que los nombres nuevos estén disponibles en CRAN ni que se hayan validado algoritmos numéricos.

## 1. Principios aprobados

1. Matemáticas antes de implementación y visualización.
2. Cuatro clases de trabajo: univariada, escalar multivariada, trayectoria y campo vectorial; la univariada es subclase del campo escalar.
3. El dominio de definición, la región de estudio, el foco visual y la ventana gráfica son objetos conceptualmente distintos.
4. Una región de **tipo III** tiene una representación de tipo I y una de tipo II del **mismo conjunto**; esta es una convención del proyecto, porque la terminología no es universal entre textos.
5. Las uniones por partes para integración aditiva se admiten cuando sus interiores son disjuntos y sus intersecciones de frontera no contribuyen a la integral considerada.
6. Las funciones existentes permanecen disponibles y con su comportamiento documentado durante la migración; cualquier incompatibilidad futura requiere aprobación y política de deprecación.
7. La discretización numérica y la visual responden a objetivos distintos; concentrar el detalle gráfico no permite desatender zonas importantes para el cálculo.
8. Plotly, ggplot2 y otros paquetes gráficos se mantienen opcionales hasta justificar un cambio mediante mediciones.
9. Los ejemplos con pipes emplean `|>` (R >= 4.1.0).

## 2. Objetos y despacho

Jerarquía de clases (de mayor a menor prioridad de despacho):

| Constructor previsto | Vector `class()` | Significado |
| --- | --- | --- |
| `vc_function1d()` | `c("vc_univariate","vc_scalar_field","vc_math")` | Aplicación real de variable real. |
| `vc_scalar()` | `c("vc_scalar_field","vc_math")` | `f: Omega subset R^n -> R`, `n >= 1`. |
| `vc_curve()` | `c("vc_trajectory","vc_math")` | `r: I -> R^m`, `m >= 1`. |
| `vc_field()` | `c("vc_vector_field","vc_math")` | `F: Omega subset R^n -> R^m`, `n,m >= 1`. |

La clase `vc_math` no implica que cualquier operador esté definido sobre todas las subclases. El método predeterminado deberá fallar de forma clara cuando una operación carezca de significado o no esté soportada.

Objetos **por composición** (no heredan de `vc_math`):

| Clase | Descripción |
| --- | --- |
| `vc_interval` | Intervalo, extremos abiertos/cerrados, orden para operaciones orientadas cuando corresponda. |
| `vc_domain` | Descripción y prueba de pertenencia al dominio, posiblemente parcial. |
| `vc_region` | Región geométrica con dimensión, representación y posibles piezas. |
| `vc_region_type1` | Region plana `a <= x <= b`, `g1(x) <= y <= g2(x)`. |
| `vc_region_type2` | Region plana `c <= y <= d`, `h1(y) <= x <= h2(y)`. |
| `vc_region_type3` | Ambas descripciones del mismo conjunto, con verificación explícita. |
| `vc_region_union` | Componentes con interiores disjuntos según contrato. |
| `vc_region_implicit` | Restricciones mediante predicados o desigualdades. |
| `vc_box` | Producto cartesiano de intervalos en dimensión `n`. |
| `vc_view` | Vista y capas gráficas, sin alterar los cálculos. |
| `vc_numerics` | Algoritmos, tolerancias y límites de trabajo numérico. |
| `vc_discretization` | Política de muestreo y refinamiento. |
| `vc_mesh` | Malla realizada: nodos, celdas, evaluaciones y diagnósticos. |
| `vc_result` | Resultado estructurado opcional con diagnósticos; no reemplaza todos los escalares por listas. |

No crear de entrada clases públicas para cada operación. Una superficie parametrizada (`vc_surface`) se contempla más adelante como objeto geométrico auxiliar para integrales de superficie, no como cuarto protagonista conceptual.

## 3. Contrato de objeto matemático

Modelo abstracto: `C = (fun, input_dim, output_dim, variables, domain, study, focus, view, numerics, metadata)`.

Campos obligatorios y invariantes:

| Campo | Tipo / invariante | Observación |
| --- | --- | --- |
| `fun` | Función R evaluable de acuerdo con un adaptador declarado | No asumir vectorización ni pureza. |
| `input_dim` | Entero >= 1 | `1` en `vc_univariate` y `vc_trajectory`. |
| `output_dim` | Entero >= 1 | `1` para todo campo escalar. |
| `variables` | Caracteres únicos, longitud `input_dim` | Evitar interpretar nombres de variables mediante evaluación de texto. |
| `domain` | Dominio de dimensión `input_dim` | Puede ser explícito, implícito o de pertenencia parcialmente verificable. |
| `study` | Región / intervalo de dimensión `input_dim` | Debe estar contenida en `domain` cuando se pueda comprobar. |
| `focus` | Región de interés visual, opcional | No redefine dominio ni región de integración. |
| `view` | `vc_view` | Solo afecta la representación. |
| `numerics` | `vc_numerics` | Configuración de métodos; no cambia la función matemática. |
| `metadata` | Lista con etiquetas, unidades opcionales, procedencia | Nunca sustituye validación. |

**Detalle importante de trayectorias:** el `domain` y `study` corresponden al intervalo del parámetro `t`, no a la imagen espacial. La orientación pertenece al intervalo parametrizado y es relevante en integrales de línea; la longitud geométrica no depende de invertirla.

**Contratos de función R:** aceptar adaptadores explícitos `input = "coordinates"` (p. ej. `function(x,y)`) o `input = "vector"` (p. ej. `function(x) x[1]+x[2]`). Detectar automáticamente solo casos inequívocos; no adivinar firmas con `...`. La prueba de evaluación verifica finitud/longitud/tipo según operación, sin requerir vectorización para todo usuario.

### 3.1 Validación

Separar:
- **Estructural:** tipos, dimensiones, cotas, integridad de listas S3, orientaciones, número de componentes.
- **Pertenencia puntual:** para `at` concreto; rechazo informativo de evaluaciones fuera de un dominio conocido.
- **Matemática:** inclusión `study subset domain`, equivalencia de las dos representaciones de tipo III, regularidad y diferenciabilidad.

Los estados internos para condiciones matemáticas son: `verified`, `sampled`, `unverified` y `invalid`. `sampled` *no constituye demostración*. El contrato debe registrar método y evidencia de cada verificación; no debe convertir silenciosamente una prueba finita en una afirmación global.

Para dominios generales descritos por funciones arbitrarias, la inclusión exacta puede ser indecidible computacionalmente. Permitir la operación con advertencia explícita cuando sea razonable, o exigir validación/hipótesis del usuario en operaciones sensibles; nunca declarar demostrada una hipótesis no verificada.

### 3.2 Comportamiento en fronteras

- Conservar extremos abiertos/cerrados y respetar singularidades.
- No evaluar automáticamente `log(x)` para `x<=0` o `1/x` en `x=0`.
- La derivada en frontera puede requerir esquema unilateral o no estar definida.
- Distinguir intervalo orientado para integrales 1D de conjunto geométrico no orientado.
- No reemplazar límites inválidos de integración por contribuciones cero silenciosamente.
- Los campos vectoriales pueden tener `input_dim != output_dim`; divergencia requiere `input_dim = output_dim` y regularidad adecuada.
- Rotacional: versión escalar en 2D, vectorial en 3D, o error/documentación cuando no haya una definición pública seleccionada.

## 4. Sistema de regiones

Convención plana:
- Tipo I: `D={(x,y):a<=x<=b, g1(x)<=y<=g2(x)}`.
- Tipo II: `D={(x,y):c<=y<=d, h1(y)<=x<=h2(y)}`.
- Tipo III: un **solo** `D` con representaciones tipo I y II válidas y equivalentes.

La convención tipo III es *local al paquete*; documentarla para evitar discrepancias con otras clasificaciones editoriales.

Constructores previstos:
- `vc_interval(lower, upper, closed = c(TRUE, TRUE))`
- `vc_rectangle(xlim, ylim)`
- `vc_disk(center, radius, boundary = "closed" | "open")`
- `vc_region_type1(xlim, lower, upper)`
- `vc_region_type2(ylim, lower, upper)`
- `vc_region_type3(by_x, by_y, validation = "unverified")`
- `vc_region_union(...)`
- `vc_region_implicit(predicate, bounds, dimension)`
- `vc_box(intervals)`

Unión por piezas: validar dimensión y relación entre componentes. Los interiores disjuntos permiten sumar integrales ordinarias cuando la intersección de las fronteras tenga medida nula y el integrando sea apropiado. No afirmar esto universalmente para integrales superficiales, objetos singulares o medidas concentradas en fronteras.

**Equivalencia tipo III:** nunca presuponer que dos representaciones dadas por funciones R son equivalentes. Documentar un mecanismo de declaración informada, certificados analíticos en familias reconocidas y verificaciones numéricas **solo exploratorias**. Una representación no verificada conserva tal estado.

Regiones implícitas y uniones no garantizan integración iterada directa: se debe elegir un algoritmo válido, descomponer la región o responder que no hay método disponible.

Para 3D planear regiones por límites iterados y cajas; no introducir una etiqueta “tipo III 3D” sin definición independiente.

## 5. Genéricos, semántica y retornos

Cada genérico recibe el objeto como primer argumento; el operador matemático no debe imprimir gráficos ni solicitar Plotly. Propuesta inicial:

| Genérico | Entrada | Retorno previsto |
| --- | --- | --- |
| `vc_evaluate(x, at, ...)` | Cualquier objeto matemático evaluable | Número o vector; múltiples puntos con estructura tabular explícita. |
| `vc_derivative(x, at = NULL, ...)` | Univariada/trayectoria, y variantes documentadas | Con `at`, derivada evaluada; sin `at`, objeto derivado cuando exista. |
| `vc_gradient(x, at = NULL, ...)` | Campo escalar | Campo vectorial derivado o valor en `at`. |
| `vc_hessian(x, at = NULL, ...)` | Campo escalar | Objeto operador / evaluación matricial según interfaz acordada. |
| `vc_jacobian(x, at = NULL, ...)` | Trayectoria/campo vectorial | Objeto operador / evaluación matricial. |
| `vc_divergence(x, at = NULL, ...)` | Campo vectorial cuadrado | Campo escalar o número evaluado. |
| `vc_curl(x, at = NULL, ...)` | Campo vectorial 2D/3D | Campo escalar o vectorial según dimensión. |
| `vc_integrate(x, ...)` | Campo escalar y región de estudio válida | Número; diagnóstico opcional mediante `details = TRUE`. |
| `vc_arc_length(x, ...)` | Trayectoria regular por partes | Número; diagnóstico opcional. |
| `vc_line_integral(field, curve, ...)` | Campo + trayectoria compatible | Número, con orientación en caso vectorial. |
| `vc_sample(x, policy = NULL, ...)` | Objeto matemático | `vc_mesh`. |
| `vc_plot(x, ...)` | Objeto / resultado / malla | Gráfico, sin efectos laterales de cálculo. |

No se promete álgebra simbólica general. Los campos derivados pueden usar evaluación numérica, siempre documentando el método y su costo. Un Hessiano en un punto devuelve matriz `n x n`; un Jacobiano de `R^n -> R^m` devuelve `m x n`.

Cuando el resultado exige diagnósticos, `details = TRUE` puede entregar `vc_result` con `value`, `method`, `estimated_error` (si existe), `evaluations`, `warnings`, `converged` y `assumptions`. Una estimación no constituye una cota certificada.

## 6. Política de discretización

**Política** `vc_discretization`, separada del objeto y de la malla; parámetros candidatos:
`strategy = c("auto","uniform","focused","adaptive")`;
`purpose = c("visualization","numerical")`;
`initial_resolution`; `max_points`; `max_depth`; `abs_tol`; `rel_tol`; `focus_weight`; `max_time` (presupuesto de trabajo, sin promesa de cronómetro exacto); `max_bytes`; `preserve` (características conocidas).

**Malla** `vc_mesh`: `nodes`, `cells`, `values`, `topology`, `purpose`, `diagnostics`, `region_ref` (referencia descriptiva, no copia costosa). Evitar exportar como compromiso público detalles internos de topología antes de estabilizarlos.

### 6.1 Tres modos

- `uniform`: rejillas/segmentos regulares; línea base y alternativa económica para problemas sencillos.
- `focused`: densidad preferente en `focus`, preservando una discretización contextual mínima.
- `adaptive`: refinamiento guiado por error local o indicador geométrico, presupuesto y condiciones de parada.
- `auto`: selecciona estrategia mediante reglas simples y benchmarks; no presuponer que adaptativo es siempre mejor.

Para visualización, una **heurística ilustrativa** de prioridad por celda es `priority(K) = eta(K) * (1 + lambda * w(K))`, con `eta` indicador normalizado y `w(K)` prioridad de foco. NO afirmar que esta prioridad representa una cota de error.

Para cálculo numérico, el refinamiento debe buscar confiabilidad integral sobre `study`; `focus` visual no puede suprimir contribuciones fuera de él. El algoritmo seleccionado debe definir su propio estimador/control de error.

### 6.2 Estrategias por geometría

- 1D: subdivisión por error de interpolación y sondeos adicionales; respetar puntos excluidos y fronteras.
- Curvas 2D/3D: refinamiento según error entre segmento y curva, geometría y parámetro; orientación y parametrización se preservan.
- Campos escalares 2D: quadtree o malla jerárquica con conectividad; comprobar que la triangulación no genere grietas.
- Campos vectoriales 2D: discretización espacial enfocada y reducción de vectores sin perder información esencial.
- Campos/sólidos 3D: postergar octree hasta justificar necesidad con benchmarks y pruebas geométricas.

**Riesgo de submuestreo:** puntos uniformes o adaptativos pueden omitir oscilaciones, picos y singularidades estrechas. Incorporar sondeos mínimos y advertencias; nunca prometer detección universal.

**Caché:** reutilizar valores compartidos dentro de la misma operación y, opcionalmente, entre operaciones únicamente con invariantes de validez. Funciones aleatorias, dependientes de estado o con efectos secundarios no deben cachearse sin consentimiento o declaración explícita de pureza. Evitar grandes cachés persistentes dentro de cada objeto.

**Rendimiento:** medir por separado tiempo de evaluación, número de evaluaciones, bytes asignados, tamaño del objeto `plotly`, tamaño serializado y tiempo de renderización cuando sea posible. Menos puntos no garantiza menos tiempo.

### 6.3 Grafismo

- 1D/2D: backend estático ligero cuando resulte apropiado, con `ggplot2` opcional.
- 3D: `plotly` opcional para interacción; no requerirlo para calcular.
- Separación `sampler -> geometry -> renderer`. Agregar trazas donde sea posible, en lugar de una por flecha.
- `vc_plot()` debe devolver el objeto gráfico; `vc_integrate()` no debe llamar `print()` ni abrir ventanas.
- Campos escalares: visualización directa cuando `n<3`. Los cálculos siguen disponibles para `n>=3`; rebanadas futuras explícitas, no proyección engañosa.

## 7. Ejemplos de API: SOLO contrato propuesto

```r
f <- vc_function1d(
  fun = function(x) sin(x) * exp(-x^2),
  domain = vc_interval(-5, 5),
  study = vc_interval(-2, 2)
)
f |> vc_derivative(at = 0)
f |> vc_integrate()
f |> vc_plot()
```

```r
Dx <- vc_region_type1(
  xlim = c(0, 1),
  lower = function(x) 0,
  upper = function(x) 1 - x
)
Dy <- vc_region_type2(
  ylim = c(0, 1),
  lower = function(y) 0,
  upper = function(y) 1 - y
)
D <- vc_region_type3(by_x = Dx, by_y = Dy)
f <- vc_scalar(
  fun = function(x, y) x + y,
  variables = c("x", "y"),
  domain = vc_rectangle(c(0, 1), c(0, 1)),
  study = D
)
# Exact reference integral: 1/3.
f |> vc_integrate(order = "dy_dx")
f |> vc_integrate(order = "dx_dy")
```

```r
curve <- vc_curve(
  fun = function(t) c(cos(t), sin(t), t / 2),
  interval = vc_interval(0, 4 * pi)
)
curve |> vc_arc_length()
curve |> vc_plot(show = c("curve", "tangent"), at = pi)
```

```r
field <- vc_field(
  fun = function(x, y) c(-y, x),
  input_dim = 2,
  output_dim = 2,
  domain = vc_rectangle(c(-2, 2), c(-2, 2))
)
field |> vc_divergence(at = c(1, 1))
```

```r
policy <- vc_discretization(
  strategy = "adaptive",
  purpose = "visualization",
  initial_resolution = 20,
  max_points = 2000
)
mesh <- f |> vc_sample(policy = policy)
```

Todos estos fragmentos son documentación de diseño, **no código ya ejecutable**.

## 8. Migración

Fases:
1. Línea base de comportamiento observable de los 47 archivos R y 46 exports declarados en `NAMESPACE`; inventario completo por símbolo, argumentos, salida y efectos secundarios.
2. Tests de regresión sobre valores, atributos, clases de retorno, warnings y gráficos cuando sea viable.
3. Nuevos constructores/métodos aislados bajo nombres `vc_*`.
4. Adaptadores internos verificables para reutilizar evaluadores y dominios.
5. Sustitución gradual de motores internos solo si se conserva comportamiento o se documenta y aprueba explícitamente el cambio.
6. Preparar una decisión de versionado después de evaluar compatibilidad.

No cambiar `DESCRIPTION`, `NAMESPACE`, `R/` ni `man/` en PR-00. No llamar "compatible" una transformación que altere sin aviso los tipos de retorno o las condiciones de integral.

## 9. Pruebas y mediciones requeridas

**Oráculos analíticos:**
- `f(x,y)=x^2+y^2`: gradiente `(2x,2y)`, Hessiano `2I`.
- `F(x,y,z)=(x,y,z)`: divergencia `3`.
- `F(x,y,z)=(-y,x,0)`: rotacional `(0,0,2)`.
- Circunferencia unitaria completa: longitud `2*pi`.
- Triángulo `x,y>=0, x+y<=1`, función `f=x+y`: integral `1/3` en ambos órdenes.

**Casos adversos:** dominios abiertos, límites singulares, funciones no vectorizadas, función no finita, orientación invertida, divisiones regionales, solapamientos no permitidos, funciones oscilatorias, discontinuidades y derivadas cerca de fronteras.

**Benchmarks obligatorios** con el mismo equipo, semilla cuando aplique y tolerancias:
1. Función univariada suave vs oscilatoria: malla uniforme vs adaptativa.
2. Superficie bivariada con foco localizado: error geométrico, evaluaciones, bytes y tiempo.
3. Curva 3D de curvatura variable: error geométrico vs puntos y costo de render.

Reportar tanto victorias como pérdidas. No asumir que `auto` supera siempre a `uniform`.

## 10. Plan de PR y criterios de aceptación

| PR | Entregable verificable |
| --- | --- |
| 00 | Arquitectura y reglas del agente; solamente documentación. |
| 01 | Tests mínimos y CI; establecer comportamiento de referencia, registrar fallos preexistentes. |
| 02 | Dominios, intervalos, regiones y validadores, con pruebas de fronteras e inclusión. |
| 03 | Constructores S3, función univariada y evaluadores, sin backends gráficos. |
| 04 | Políticas de discretización y muestreo 1D, con benchmarks. |
| 05 | Muestreo adaptativo 2D y de curvas, conectividad y presupuestos verificados. |
| 06 | Operadores diferenciales y resultados con pruebas analíticas. |
| 07 | Integración y campos vectoriales, regiones tipo III y uniones por partes. |
| 08 | Backends gráficos desacoplados con mediciones de tamaño, velocidad y memoria. |
| 09 | Adaptadores y tests de regresión para API heredada, sin cambios silenciosos. |
| 10 | Contenido pedagógico, documentación formal, pkgdown y verificaciones CRAN. |

Cada PR debe incluir pruebas nuevas relevantes, documentación de parámetros, ejemplos reproducibles, justificación matemática y observaciones de compatibilidad. No fusionar por cobertura numérica agregada únicamente.

## 11. Preguntas que requieren decisión antes de codificar

1. `vc_gradient(x)` devuelve campo derivado, y `vc_gradient(x, at=...)` valor: confirmar homogeneidad con `vc_hessian`/`vc_jacobian`, cuyos objetos operacionales pueden requerir un tipo auxiliar.
2. ¿El dominio desconocido queda como `unverified` por defecto o se exige al usuario declarar una región que actúe como dominio de trabajo?
3. ¿Cómo declara el usuario la equivalencia tipo III sin fingir verificación automática? Propuesta: etiqueta de estado y evidencia opcional.
4. ¿Qué semántica exacta tiene una operación con punto válido en dominio pero fuera de `study`? Propuesta: evaluación permitida sin alterar el estudio.
5. ¿Qué resultados de integración retornan por defecto valores simples y cuáles incluyen diagnóstico?
6. ¿Qué características se conservan como efectos secundarios en las funciones heredadas que actualmente imprimen gráficos?
7. ¿Usar nombres largos `vc_function1d`/`vc_region_type1` para API estable y reservar alias cortos solo tras pruebas docentes?

Estas preguntas no invalidan el diseño. Deben resolverse por revisión explícita, no por decisiones silenciosas de Codex.

## 12. Referencias técnicas

- Wickham, *Advanced R*, capítulo S3: https://adv-r.hadley.nz/s3.html
- R Core, *Writing R Extensions*, registro de métodos S3 y namespaces: https://cran.r-project.org/doc/manuals/r-release/R-exts.html
- R Core, `stats::integrate`: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/integrate.html

Estas fuentes respaldan las convenciones de R; la arquitectura numérica de este documento es una **propuesta del proyecto**, no algo implementado o avalado por las referencias.
