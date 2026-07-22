# Criterios de inclusion de tripletes

Guia breve de los requisitos de evidencia que debe cumplir un triplete antes de
cargarlo en `input/ddi_reference_input.xlsx`. El set sirve como **set de
referencia** de interacciones droga-droga (IDD) pediatricas para el benchmark, y
prioriza especificidad y trazabilidad.

Cada triplete es un **control positivo** o un **control negativo** (`control_type`).
El dataset está armado en linea con Kontsioti et al.,
*Pharmacoepidemiol Drug Saf* 2023;32:832-844 (que muestra que el marco AUC/PPV
exige controles positivos y negativos, y que toda restriccion del set de
evaluacion debe justificarse explicitamente).

## Que es un triplete

Una unidad `farmaco1 x farmaco2 x evento MedDRA`: dos farmacos coadministrados y
un evento adverso. En un **control positivo** el evento es atribuible a la
**interaccion** del par (no a un solo farmaco). En un **control negativo** el par
se coadministra de forma plausible pero **no hay interaccion documentada** que
vincule al par con ese evento. Un mismo par puede generar varios tripletes, uno
por evento.

## Requisitos minimos (todos obligatorios)

1. **Poblacion pediatrica.** La evidencia debe provenir de pacientes < 21 años
   (neonato a adolescente). No se aceptan extrapolaciones desde adultos salvo que
   exista confirmacion pediatrica directa del par.
2. **Par especifico.** La fuente debe documentar la coadministracion de **ambos**
   farmacos, no cada uno por separado. Ambos deben resolver a un codigo ATC 5th
   del vocabulario (se eligen del desplegable).
3. **Evento atribuible a la interaccion.** Debe haber plausibilidad temporal y/o
   mecanistica de que el evento surge de la interaccion (p.ej. toxicidad por
   sobreexposicion, perdida de eficacia por induccion, precipitacion fisico-quimica).
4. **Mecanismo declarado.** `interaction_type` debe ser uno de:
   `pharmacokinetic`, `pharmacodynamic`, `mixed`, `pharmaceutical` o `unknown`.
   Usar `unknown` solo si el mecanismo no esta resuelto en la literatura. Los
   controles negativos usan `none` (no tienen mecanismo de interaccion).
5. **Evento mapeable a MedDRA.** El evento se carga en su nivel mas fino
   disponible (LLT/PT/HLT/HLGT) eligiendolo del desplegable. Si el termino exacto
   no existe como PT, mapear al PT mas cercano y dejar constancia en `comments`.
6. **Al menos una fuente trazable.** Cada `triplet_id` debe tener >= 1 fila en la
   hoja `sources` con PMID o DOI (o, en su defecto, URL estable) y la cita
   completa.
7. **Nivel de evidencia declarado.** `evidence_level` debe ser uno de los niveles
   controlados (ver mas abajo). Si se deja vacio, el script lo deriva del nivel
   mas alto de las fuentes; el curador puede sobreescribirlo.
8. **Detectabilidad en FAERS.** El triplete completo
   `farmaco1 x farmaco2 x evento` debe tener **>= 1 coReporte pediatrico** en
   FAERS (el par y el evento coocurren en al menos un caso). Sin coReporte del
   triplete la celda co-expuesta-con-evento del benchmark queda vacia: un positivo
   nunca puede ser verdadero positivo y un negativo es un verdadero negativo
   trivial que infla la especificidad. El coReporte del *par* (req. 2) es
   necesario pero **no** suficiente: hay que verificar el evento curado. Chequear
   con `scripts/R/faers_triplet_coreport.R` (imprime el total y el desglose por
   etapa NICHD sin cargar el dataset) y rechazar el triplete si el total es 0. El
   conteo se hace a nivel **MedDRA PT** (el mas fino): es el gate estricto, porque
   un coReporte a PT tambien cuenta en su HLT/HLGT, de modo que >= 1 a PT garantiza
   detectabilidad en cualquier nivel de roll-up en que corra el benchmark.

## Niveles de confianza (`confidence_level`)

Solo se admiten dos valores; elegir el mayor que la evidencia sustente:

- **high**: etiqueta regulatoria (FDA/EMA), ensayo controlado, meta-analisis o
  revision sistematica; **o** caso/serie pediatrica con relacion temporal y
  mecanistica clara, idealmente con confirmacion independiente.
- **moderate**: caso unico, evidencia de nivel resumen (revision narrativa),
  mecanismo incierto, o el evento es plausible pero la atribucion es menos
  directa.

No se incluyen tripletes que solo alcanzarian un nivel "bajo" (interaccion
meramente teorica, sin caso pediatrico ni respaldo regulatorio).

## Nivel de evidencia (`evidence_level`)

Campo controlado que clasifica la **fuerza de la evidencia** que sustenta la
interaccion (o, en un negativo, la fuente que documenta su **ausencia**). Replica
la taxonomia BNF/Micromedex usada en la literatura de sets de referencia de IDD,
para poder examinar de forma reproducible el impacto del nivel de evidencia sobre
la evaluacion (Kontsioti 2023). De mayor a menor:

- `regulatory_label`: la etiqueta FDA/EMA o el SmPC reconoce (o, en un negativo,
  omite) la interaccion.
- `controlled_study_or_meta`: ECA, meta-analisis, revision sistematica o estudio
  poblacional/PK de modelado.
- `observational_study`: estudio retrospectivo, clinico o de cohorte/PK.
- `case_series`: serie de casos.
- `single_case_report`: reporte de caso unico.
- `theoretical`: solo mecanismo o extrapolacion de clase, sin caso observado
  (excluido por los criterios de los positivos; se mantiene por completitud).

Si el curador deja `evidence_level` vacio, la consolidacion lo deriva tomando el
nivel mas alto presente entre el `source_type` del triplete y todas sus fuentes
(`derive_triplet_evidence_level` en `00_functions.R`). Es un default determinista
y reproducible que el curador puede sobreescribir.

## Origen de los candidatos positivos

Los positivos parten de una lista externa para que la seleccion del par no quede sesgada por el agente.
`scripts/R/01_generate_positive_candidates.R` interseca los controles positivos de
CRESCENDDI (Kontsioti et al., *Sci Data* 2022;9:72; tiers Micromedex
`Established`/`Probable`) con el coReporte pediatrico FAERS, y emite una lista
ordenada por cobertura de etapas NICHD, `pair_coreport` y nivel de evidencia.
CRESCENDDI ancla el *par* (interaccion documentada en adultos) y aporta el evento
adulto **solo como pista**; el evento pediatrico, su atribucion a la interaccion y
el mapeo MedDRA se curan igual que cualquier positivo (requisitos de arriba).

Cobertura FAERS por etapa: el par debe estar coadministrado en **>= 2 etapas
NICHD**. No se exige cobertura total. Esta cobertura es del *par*; una vez elegido
el evento pediatrico, verificar ademas la detectabilidad del **triplete** (req. 8)
con `scripts/R/faers_triplet_coreport.R`, porque `01` cuenta el par, no el evento.

## Controles negativos

Un control negativo (`control_type = negative`) es un triplete
`farmaco1 x farmaco2 x evento` que sirve como **verdadero negativo** para estimar
especificidad. Requisitos:

1. **Par plausible y triplete detectable.** Ambos farmacos deben resolver a un
   ATC 5th y ser coadministrables en pediatria, y el triplete completo debe tener
   **>= 1 coReporte pediatrico** en FAERS (par + evento en el mismo caso; req. 8),
   para que el negativo sea informativo y no un verdadero negativo trivial.
   `scripts/R/03_generate_negative_candidates.R` ya impone este piso
   (`triplet_coreport >= 1`, columna `triplet_coreport` en la planilla).
2. **Sin interaccion documentada.** No debe existir interaccion conocida del par
   para ese evento en los compendios/etiquetas consultados. `interaction_type` es
   `none` y `mechanism` queda vacio.
3. **Ausencia trazable.** Igual que un positivo, requiere >= 1 fila en `sources`,
   pero la cita documenta la **ausencia de evidencia** (p.ej. "sin interaccion
   listada en DDINTER 2.0/FDA/SmPC del par al <fecha>"), con la fecha de consulta.
4. **Evitar misclasificacion.** No usar como negativo un par cuyo evento sea
   atribuible a un solo farmaco, ni pares para los que solo falta evidencia por
   ser poco estudiados (ausencia de evidencia != evidencia de ausencia). Esta
   cautela responde a la misclasificacion de negativos documentada en la
   literatura (Hauben 2016).
5. **Sin modulacion ontogenica.** `ontogenic_modulation` es `no` por definición.

`confidence_level` en un negativo expresa la confianza en la **clasificacion como
no-interaccion** (misma escala `high`/`moderate`).

Para no construir los negativos a mano desde cero, `scripts/R/03_generate_negative_candidates.R`
recombina el set positivo con dos estrategias apareadas (`event_swap` y `drug_swap`,
ambas presentes) y usa el coReporte pediatrico real en FAERS (del par y del
triplete) **solo como piso de elegibilidad**, nunca para rankear. Entre los
candidatos elegibles la seleccion es **aleatoria con semilla fija**
(`negative_seed`), reproducible y sin sesgo por cantidades derivadas del desenlace
(Kontsioti 2022: los negativos se sortean sin agregar sesgos por diseño). No hay
columna `suggested`.

`single_drug_event_max` (veces que el evento coocurre con un solo farmaco) se
reporta **solo para estratificar** a posteriori (Kontsioti sugiere estratificar por
el ADR mono-farmaco, no filtrar por el). **No** es criterio de descarte: la
atribucion mono-farmaco (paso 4) se evalua contra etiquetas/compendios, no
umbralando ese conteo de FAERS, que comparte celdas con el AC y su uso como filtro
seria circular.

Esos candidatos automatizan el paso 1 (par plausible y detectable), pero **no** son
negativos hasta que el curador verifica la ausencia documentada (pasos 2-3) y la no
atribucion mono-farmaco (paso 4) contra las fuentes, y los carga en el workbook.

## Modulacion ontogenica (opcional)

Si la evidencia indica que el riesgo de la IDD **varia con la etapa del
desarrollo**, registrarlo para poder contrastarlo con la salida por etapa del
benchmark:

- `ontogenic_modulation`: `yes` si hay evidencia de modulacion por edad, `no` si
  la evidencia indica que no la hay, `unknown` (default) si no se evaluo.
- `higher_risk_stages`: etapa(s) NICHD donde el riesgo es mayor, como subconjunto
  separado por comas de:
  `term_neonatal, infancy, toddler, early_childhood, middle_childhood,
  early_adolescence, late_adolescence`.
  Solo se completa cuando `ontogenic_modulation = yes`.
- `ontogeny_evidence`: breve justificacion del patron etario (mecanismo,
  ventana de riesgo).

Estas etapas deben coincidir con `niveles_nichd` (definido en `00_functions.R`),
que a su vez replica el orden usado por `gam_benchmark`.

## Exclusiones

Aplican a los **controles positivos** (los negativos se rigen por su seccion):

- Evidencia solo en adultos, sin confirmacion pediatrica del par.
- Interacciones puramente teoricas o predichas por bases de datos sin caso real.
- Evidencia únicamente por detección de señales de desproporcionalidad (evitar circularidad en el benchmark).
- Eventos atribuibles a un unico farmaco (no a la interaccion).
- Eventos que no pueden mapearse a ningun nivel MedDRA del vocabulario local.
- Farmacos sin codigo ATC 5th en el vocabulario.

## Como agregar un triplete

1. Abrir `input/ddi_reference_input.xlsx`, hoja `triplets`, y completar una fila
   nueva usando los desplegables. Elegir `control_type` (`positive`/`negative`),
   los farmacos, el evento, `interaction_type` (`none` si es negativo),
   `evidence_level`, `confidence_level` y las etapas.
2. Agregar su(s) fuente(s) en la hoja `sources` con el mismo `triplet_id` (en un
   negativo, la cita documenta la ausencia de interaccion y la fecha de consulta).
3. Correr `scripts/R/curate_pediatric_ddi_reference_set.R`. El script valida
   el mapeo y las restricciones; si algo no cumple, falla indicando el problema.

> Si el workbook se regenera desde el seed (`00_build_input_template.R` con
> `overwrite_existing <- TRUE`), las nuevas columnas (`control_type`,
> `evidence_level`) ya vienen con sus desplegables. Sobre un workbook viejo,
> `01` completa `control_type = positive` y deriva `evidence_level`
> automaticamente, de modo que el pipeline sigue corriendo sin regenerar.
