# Flujo de trabajo de curacion (SOP)

Procedimiento operativo para curar el set de referencia pediatrico de IDD,
**tanto controles positivos como negativos**, de forma reproducible y trazable.
Criterios de aceptacion en
[`INCLUSION_CRITERIA.md`](INCLUSION_CRITERIA.md) y la estructura en
[`README.md`](README.md). Seguir este orden evita reprocesos y mantiene cada
triplete auditable.

## Principio rector

La **unica entrada manual** es el workbook `input/ddi_reference_input.xlsx`: ahi
el curador humano aprueba cada triplete que entra al set. Todo lo demas (mapeo a
vocabulario, validacion, entregables de `results/`, candidatos positivos y
negativos) es **generado por script y reproducible**: borrar `results/` y volver
a correr los scripts reconstruye el set identico desde el workbook. No se edita
codigo ni `results/` para agregar tripletes, ni se transcriben numeros a mano.

## Mapa del flujo

```
                 (manual)     (curate_pediatric_ddi_reference_set)
  vocabulario --> workbook  ----------------------> curated_*.csv --> gam_benchmark
   (dropdowns)   triplets +                          (entregable)
                 sources

  POSITIVO: candidatos (script 01) --> agente: dossier ---+
            [CRESCENDDI + FAERS]      + filas workbook     |--> fila en workbook
  NEGATIVO: candidatos (script 03) --> tamiz manual -------+     + fila en sources
            [recombinacion + FAERS]   [compendios]
```

`curate_pediatric_ddi_reference_set` es el paso de **consolidacion**: corre
despues de cada edicion del workbook y es lo unico que escribe los entregables
de `results/`. El script `03` parte de los positivos ya consolidados, asi que la
consolidacion aparece dos veces en el ciclo completo:

```
  01_generate_positive_candidates
    --> curacion manual de positivos (seccion 1)  --> workbook
      --> curate_pediatric_ddi_reference_set      (seccion 3)
        --> 03_generate_negative_candidates
          --> tamiz manual de negativos (seccion 2) --> workbook
            --> curate_pediatric_ddi_reference_set  (seccion 3)
              --> gam_benchmark
```

## 0. Setup (una sola vez)

Construir el workbook con sus desplegables (no sobrescribe uno existente; protege
la curacion manual):

```powershell
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' scripts\R\00_build_input_template.R
```

El agente de curacion espera el conector MCP `biomcp` (declarado en `.mcp.json`),
que aporta etiquetas FDA (`openfda_label_searcher`/`openfda_label_getter`), FAERS
(`openfda_adverse_searcher`), PubMed (`article_searcher`) y ensayos
(`trial_searcher`); PubMed tambien via el MCP `search_articles`. La DB de
interacciones es **DDInter local**: `scripts/R/ddinter_lookup.R` sobre los CSV en
`input/ddinter/` (pair-level; ausencia del par = evidencia de negativo). El
vocabulario del workbook se consulta con `scripts/R/vocab_lookup.R` (devuelve solo
coincidencias, sin cargar las hojas de referencia enteras).

## 1. Curar un control POSITIVO

La curacion de positivos **parte de una lista de candidatos**, no de pares
elegidos libremente, para que la seleccion del par sea externa y reproducible (no
sesgada por el modelo).

### 1a. Generar candidatos (automatico)

```powershell
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' scripts\R\01_generate_positive_candidates.R
```

Produce `results/positive_control_candidates/positive_control_candidates.csv`:
pares de CRESCENDDI (tiers Micromedex `Established`/`Probable`) coadministrados en
FAERS pediatrico, ordenados por cobertura de etapas NICHD, `pair_coreport` y nivel
de evidencia. Columnas de apoyo: `crescenddi_event_hint` (evento adulto, solo
pista), `coadmin_<etapa>`, `stages_covered`, `meets_coverage`. El requisito de
cobertura es **>= 2 etapas NICHD** (`term_neonatal` esta estructuralmente vacio
para los pares de CRESCENDDI). Requiere la planilla CRESCENDDI ya descargada en
`dual_curation_validation/input/` (ver ese README).

### 1b. Curar cada candidato (agente -> humano)

Trabajar la lista de arriba hacia abajo hasta alcanzar el objetivo (50
positivos). Por cada candidato, el **agente**: busca la evidencia pediatrica
(skill `ddi-positive-curation`), escribe el dossier en `suggested/positivos/`,
verifica la coadministracion FAERS por etapa y redacta el borrador de las filas
del workbook. El **humano cura** ese borrador (acepta/edita/descarta).

1. **Verificar el umbral de evidencia** contra
   [`INCLUSION_CRITERIA.md`](INCLUSION_CRITERIA.md) (poblacion pediatrica, par
   especifico coadministrado, evento atribuible a la *interaccion*, mecanismo,
   evento mapeable, >=1 fuente trazable). CRESCENDDI ancla el *par*; el evento
   pediatrico y su atribucion se curan (no se hereda el evento adulto). Si no
   cumple, no se carga.
2. **Hoja `triplets`** -> nueva fila usando los desplegables:
   - `control_type = positive`.
   - `drug1` / `drug2` desde el desplegable ATC 5th (`sustancia; via`).
   - evento en su nivel MedDRA mas fino disponible (`event_llt`/`event_pt`/
     `event_hlt`/`event_hlgt`); el script rellena los niveles mas gruesos.
   - `interaction_type` (`pharmacokinetic`/`pharmacodynamic`/`mixed`/
     `pharmaceutical`/`unknown`, nunca `none`), `mechanism`, `confidence_level`
     (`high`/`moderate`). `evidence_level` opcional (si se deja vacio, el script
     deriva el nivel mas alto de las fuentes).
   - `ontogenic_modulation` y, si es `yes`, `higher_risk_stages` (etapas NICHD).
3. **Hoja `sources`** -> >=1 fila con el mismo `triplet_id` (PMID/DOI o URL
   estable + cita completa).
4. **Correr la consolidacion** (seccion 3) para mapear y validar.

## 2. Curar un control NEGATIVO

Un negativo es un par plausible *sin interaccion documentada* para ese evento
(ver [`INCLUSION_CRITERIA.md` -> Controles negativos](INCLUSION_CRITERIA.md)). Se busca evitar errores de **misclasificacion** ("ausencia de evidencia != evidencia de
ausencia", Hauben 2016)..

### 2a. Generar candidatos (automatico)

```powershell
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' scripts\R\03_generate_negative_candidates.R
```

Produce `results/negative_control_candidates/negative_control_candidates.csv`:
recombinaciones apareadas del set positivo con dos estrategias (`event_swap` y
`drug_swap`, ambas presentes). El coReporte pediatrico real en FAERS (del par y del
triplete) se usa **solo como piso de elegibilidad**; entre los elegibles la lista
sale **ordenada al azar con semilla fija** (`negative_seed`), sin columna
`suggested`. Trabajar la lista de arriba hacia abajo, sin reordenar: ese orden es la
seleccion reproducible y sin sesgo de Kontsioti. Columnas a leer por fila (contexto,
no criterios de descarte):

- `match_strategy` + `known_interacting_pair`: `drug_swap` = par nuevo (`FALSE`);
  `event_swap` con `TRUE` = el par interactua para *otro* evento, confirmar que NO
  interactua para *este*.
- `pair_coreport`: coReportes pediatricos del par (plausibilidad de coadministracion).
- `single_drug_event_max`: veces que el evento coocurre con UN solo farmaco. **Solo
  una pista** para revisar la etiqueta (y variable de estratificacion posterior),
  **no** un umbral de descarte: la atribucion mono-farmaco se decide por la etiqueta.
- `triplet_coreport >= 1`: piso de detectabilidad ya aplicado por el script.

### 2b. Tamiz manual 

**Descartar si**:

1. Existe interaccion del par documentada para ese evento en compendios/etiquetas
   (consultar **FDA + SmPC + una DB de IDD** como DDInter/Drugs.com).
2. El evento es un ADR conocido de *cualquiera* de los dos farmacos por separado
   (decidir por la **etiqueta** de cada farmaco; `single_drug_event_max` es solo la
   pista de FAERS que motiva la consulta, no la prueba).
3. El par casi no se coadministra en pediatria o solo "falta evidencia" por ser
   poco estudiado (no es verdadero negativo).

### 2c. Cargar el negativo aceptado

1. **Hoja `triplets`**: copiar las primeras 23 columnas de la fila del candidato
   (ya vienen con los valores de desplegable correctos). Asignarle un
   **`triplet_id` estable y permanente** propio (p.ej. `N001`, `N002`, ...) — **no
   reutilizar el id `Nxxx` del CSV**, que se renumera en cada corrida (ver
   seccion 4). Dejar `control_type = negative`, `interaction_type = none`,
   `mechanism` vacio, `ontogenic_modulation` distinto de `yes`, `confidence_level`
   (`high`/`moderate`, confianza en la clasificacion como no-interaccion).
2. **Hoja `sources`**: >=1 fila con ese `triplet_id` que documenta la **ausencia**
   y la **fecha de consulta** (p.ej. cita: "sin interaccion listada en
   BNF/Micromedex/SmPC del par al AAAA-MM-DD"), con URL estable.
3. **Correr la consolidacion** (seccion 3).

## 3. Validar y consolidar

```powershell
# Mapea, valida restricciones y escribe los entregables curados.
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' scripts\R\curate_pediatric_ddi_reference_set.R
```

El script falla con un mensaje claro si algo no cumple. Causas tipicas:

- *farmaco no coincide con ATC 5th* -> revisar que el valor venga del desplegable.
- *negativo con `interaction_type` != `none`* / *positivo con `none`* -> corregir.
- *negativo con `ontogenic_modulation = yes`* -> un negativo no modula la IDD.
- *triplete sin fuente* / *triplete semantico duplicado* -> completar `sources` o
  eliminar el duplicado.

Cuando pasa, regenerar el benchmark downstream:

```powershell
cd ..\gam_benchmark
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' scripts\R\01_validate_pediatric_ddi_benchmark.R
```

Con negativos cargados, las metricas dependientes de negativos (especificidad,
PPV, NPV, AUC) se pueblan automaticamente; con set positivo-only quedan `NA`.

## 4. Reproducibilidad y documentacion

- **`triplet_id` estable y unico**: una vez asignado no se reusa ni se reordena,
  aunque se elimine el triplete. Los ids `Nxxx` del CSV de candidatos son
  **posicionales** (cambian entre corridas): al aceptar, fijar un id propio.
- **Toda fila trazable**: positivos con PMID/DOI/URL de la evidencia; negativos
  con la cita de ausencia + fecha de consulta. Sin fuente, la consolidacion falla.
- **No transcribir numeros**: los outputs y el manuscrito leen los CSV en vivo.
- **Reproducibilidad verificable**: borrar `results/` y correr la consolidacion
  (y `01`, `03`) reconstruye todo desde el workbook; `01` es determinista y `03`
  usa una semilla fija (`negative_seed`), asi que la lista de candidatos negativos
  es identica en cada corrida.


## Checklist 

**Positivo:** correr `01` -> tomar candidato de la lista (cobertura >= 2 etapas)
-> evidencia pediatrica + dossier -> criterio OK -> fila `triplets`
(`control_type=positive`, `interaction_type` real) -> fila(s) `sources`
(PMID/DOI) -> humano cura -> correr la consolidacion.

**Negativo:** correr `03` -> tomar candidatos en el orden aleatorio del CSV ->
tamiz compendios (sin interaccion para ese evento + no es ADR mono-farmaco por
etiqueta) -> fila `triplets` (`control_type=negative`, `interaction_type=none`, id
estable propio) -> fila `sources` (ausencia + fecha) -> correr la consolidacion.
