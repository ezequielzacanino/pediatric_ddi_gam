# DDI Reference Set

Construye el set de referencia curado pediatrico de interacciones droga-droga. Cada
triplete es un **control positivo** (interaccion curada) o **negativo** (par/evento
sin interaccion documentada), con evidencia específica en población pediátrica. Produce el dataset curado que
consume `gam_benchmark/`. Los criterios de evidencia (incluidos los del nivel de
evidencia y los controles negativos) estan en
[`INCLUSION_CRITERIA.md`](INCLUSION_CRITERIA.md); el procedimiento paso a paso para
curar positivos y negativos esta en [`CURATION_WORKFLOW.md`](CURATION_WORKFLOW.md).

## Orden de ejecucion

Desde la raiz de `ddi_reference_set/`:

```powershell
# Construye (o regenera) workbook .xlsx vacio con sus dropdowns (sin precargar tripletes).
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' scripts\R\00_build_input_template.R

# Genera plantilla de candidatos positivos a partir de CRESCENDDI intersecados con coReporte pediatrico FAERS.
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' scripts\R\01_generate_positive_candidates.R

# Consolida plantilla. Mapea y produce los CSV de salida. Corre despues de curar positivos y otra vez despues de curar negativos.
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' scripts\R\curate_pediatric_ddi_reference_set.R

# Genera plantilla de candidatos negativos, a partir de recombinación del set positivo y chequeo de presencia en FAERS.
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' scripts\R\03_generate_negative_candidates.R
```

El orden de corrida completo esta en
[`CURATION_WORKFLOW.md`](CURATION_WORKFLOW.md) (seccion "Mapa de flujo").

## Como agregar tripletes

La curacion agéntica se hace en `input/ddi_reference_input.xlsx` (hoja `triplets`). Cada fila es un triplete:

- `control_type`: `positive` (interaccion curada) o `negative` (par/evento sin
  interaccion documentada). En un negativo, `interaction_type = none`.
- `evidence_level`: nivel de evidencia controlado (ver `INCLUSION_CRITERIA.md`).
  Si se deja vacio, el script lo deriva del nivel mas alto de las fuentes.
- `drug1` / `drug2`: se eligen del desplegable de nombres ATC 5th (formato
  `sustancia; via`), conectado al vocabulario para evitar errores de tipeo.
- `event_llt` / `event_pt` / `event_hlt` / `event_hlgt`: el evento adverso se
  carga en el nivel MedDRA mas fino disponible, cada columna con su propio
  desplegable. El script toma el nivel mas fino completado y rellena
  automaticamente los niveles mas gruesos (PT/HLT/HLGT). Las columnas mas finas
  que el nivel ingresado quedan vacias.
- `interaction_type` (incluye `pharmaceutical`) y `confidence_level` tienen desplegable.
- `ontogenic_modulation` / `higher_risk_stages` / `ontogeny_evidence`: registran si
  el riesgo de la IDD esta modulado por la edad y en que etapa(s) NICHD es mayor,
  para poder contrastarlo con la deteccion por etapa del benchmark.
  `higher_risk_stages` es un subconjunto separado por comas de las 7 etapas
  canonicas (`niveles_nichd`, alineadas con `gam_benchmark`) y solo se completa
  cuando `ontogenic_modulation = yes`.

Los criterios de evidencia para aceptar un triplete estan en
[`INCLUSION_CRITERIA.md`](INCLUSION_CRITERIA.md).

Las fuentes se cargan en la hoja `sources`, una fila por cita, referidas por
`triplet_id`. El script `00` no sobrescribe la planilla existente (protege la
curacion manual); para regenerarla vacia hay que poner
`overwrite_existing <- TRUE`.

El dataset final que se usará en `gam_benchmark` es el de `results/curated_pediatric_ddi_reference_set/`, el cual solo debe ser editado por curador humano. 

## Estructura

- `00_functions.R`: helpers de mapeo MedDRA. `build_vocabulary_picklists()` arma
  los desplegables y `resolve_meddra_event_levels()` resuelve el evento al nivel
  mas fino ingresado y lo consolida a PT/HLT/HLGT (LLT se traduce a su PT via la
  relacion "Is a").
- `scripts/R/00_build_input_template.R`: genera `input/ddi_reference_input.xlsx`
  vacio, con solo la fila de encabezado y sus desplegables.
- `scripts/R/01_generate_positive_candidates.R`: interseca los controles positivos
  de CRESCENDDI (tiers Micromedex `Established`/`Probable`) con el coReporte
  pediatrico FAERS y emite `results/positive_control_candidates/` con la cobertura
  por etapa NICHD y el evento adulto como pista. Es la lista-guia de la que parte
  la curacion de positivos (en vez de elegir pares de forma libre).
- `scripts/R/curate_pediatric_ddi_reference_set.R`: lee la planilla, resuelve
  ATC y eventos contra el vocabulario y escribe los outputs. Es el paso de
  consolidacion: corre despues de cada edicion del workbook, tanto al curar
  positivos como negativos.
- `scripts/R/03_generate_negative_candidates.R`: recombina el set
  positivo (matched 1:1: `event_swap`/`drug_swap`), filtra por coReporte
  pediatrico real en FAERS y emite `results/negative_control_candidates/` con
  evidencia de plausibilidad y flags de atribucion mono-farmaco. Automatiza los
  pasos 1-2 del proceso; el tamiz contra compendios y la aceptacion siguen siendo
  curacion manual.
- `agent/`: carpeta autocontenida del agente de curacion (contexto + trabajo);
  ver [`agent/README.md`](agent/README.md).
  - `agent/skills/`: procedimiento paso a paso por tipo de control
    (`SKILL.md` + `TEMPLATE.md` de positivos y negativos).
  - `agent/tools/`: helpers que corre el agente durante la curacion
    (`vocab_lookup.R`, `ddinter_lookup.R`, `faers_triplet_coreport.R`).
  - `agent/workspace/{positivos,negativos}/`: dossiers de evidencia que produce
    el agente (local-only).
- `input/`: planilla de entrada curada a mano.
- `results/`: outputs generados.

## Entradas esperadas

El vocabulario OMOP vive en una raiz compartida del workspace, referenciada por
ruta relativa:

- `../data/vocabulary/vocabulary_SNOMED_MEDDRA_RxNorm_ATC/CONCEPT.csv`
- `../data/vocabulary/vocabulary_SNOMED_MEDDRA_RxNorm_ATC/CONCEPT_ANCESTOR.csv`
- `../data/vocabulary/vocabulary_SNOMED_MEDDRA_RxNorm_ATC/CONCEPT_RELATIONSHIP.csv`
  (solo se lee cuando la planilla tiene eventos cargados a nivel LLT)

## Outputs (entregables a gam_benchmark)

- `results/curated_pediatric_ddi_reference_set/curated_pediatric_ddi_triplets.csv`
- `results/curated_pediatric_ddi_reference_set/curated_pediatric_ddi_sources.csv`

`gam_benchmark/` los lee por ruta relativa
(`../ddi_reference_set/results/curated_pediatric_ddi_reference_set/`).
