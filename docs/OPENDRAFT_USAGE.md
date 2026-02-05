# OpenDraft (Engine) - How to Run and Use

This document explains how to run the OpenDraft engine that lives in:
`/home/dodo/Documents/projects/opendraft/engine`

It covers environment activation, required `.env.local` settings, and the main commands from fastest to slowest.

## 1) Activate the environment
Use one of the following:

```bash
cd /home/dodo/Documents/projects/opendraft/engine
source .venv/bin/activate
```

or

```bash
conda activate pharmacophore-toolkit
cd /home/dodo/Documents/projects/opendraft/engine
```

## 2) One-time configuration (.env.local)
Create/modify: `/home/dodo/Documents/projects/opendraft/engine/.env.local`

```bash
MODEL_PROVIDER=openai
MODEL_NAME=deepseek-v3.2-thinking
OPENAI_BASE_URL=https://llm.ai.e-infra.cz/v1
OPENAI_API_KEY=YOUR_EINFRA_TOKEN

# Optional Google Gemini (for grounded search)
GOOGLE_API_KEY=YOUR_GOOGLE_AI_STUDIO_KEY

# Optional: reduce 429s from Semantic Scholar
SEMANTIC_SCHOLAR_API_KEY=YOUR_S2_KEY

# Optional: PubMed etiquette
NCBI_EMAIL=you@example.com
NCBI_TOOL=opendraft

# Optional: source toggles
ENABLE_PUBMED=true
ENABLE_ARXIV=true
ENABLE_GEMINI_GROUNDED=true

# Optional: rate limits (safe defaults)
SCOUT_PARALLEL_WORKERS=1
SCOUT_BATCH_SIZE=1
SCOUT_BATCH_DELAY=4
RATE_LIMIT_DELAY=2

# Optional: deep research cap
MAX_RESEARCH_QUERIES=25

# Optional: PDF download + summary
DOWNLOAD_PDFS=false
SUMMARIZE_PDFS=false
MAX_PDF_DOWNLOADS=10
PDF_SUMMARY_MAX_CHARS=12000

# Optional: OpenAI-compatible timeouts/retries
OPENAI_TIMEOUT=240
OPENAI_CONNECT_TIMEOUT=10
OPENAI_RETRIES=3
OPENAI_RETRY_BACKOFF=2.0
MODEL_MAX_OUTPUT_TOKENS=1200
```

## 3) Commands (fastest to slowest)

### A) Fastest: LLM smoke test
Checks API connectivity and model config.

```bash
cd /home/dodo/Documents/projects/opendraft/engine
python - <<'PY'
from utils.agent_runner import setup_model
m = setup_model()
print(m.generate_content("Say hello in one sentence."))
PY
```

### B) Fast: search-only, no deep research
Scout + summaries, skips deep research planning.

```bash
python draft_generator.py \
  --topic "3D pharmacophore modeling in ligand based drug design" \
  --academic-level research_paper \
  --search-only \
  --no-deep-research \
  --skip-citation-gate
```

### C) Medium: search-only with capped deep research
Deep research planning but limited query count.

```bash
python draft_generator.py \
  --topic "3D pharmacophore modeling in ligand based drug design" \
  --academic-level research_paper \
  --search-only \
  --max-research-queries 25 \
  --skip-citation-gate
```

### D) Slow: search + download PDFs + summarize
Best-effort OA download, summaries appended to research output.

```bash
python draft_generator.py \
  --topic "3D pharmacophore modeling in ligand based drug design" \
  --academic-level research_paper \
  --search-only \
  --max-research-queries 25 \
  --skip-citation-gate \
  --download-pdfs \
  --summarize-pdfs
```

### E) Slowest: full draft
Runs the full pipeline and writes draft artifacts.

```bash
python draft_generator.py \
  --topic "3D pharmacophore modeling in ligand based drug design" \
  --academic-level research_paper
```

## 4) Where outputs go
All outputs are saved under:

```
/home/dodo/Documents/projects/opendraft/engine/tests/outputs/generated_draft/
```

Key files:
- `research/scout_raw.md` — raw citation results
- `research/combined_research.md` — summarized research
- `research/research_gaps.md` — identified gaps
- `research/papers/` — per‑paper summaries
- `research/pdfs/` — downloaded PDFs + summaries (if enabled)
- `drafts/` — outlines and draft sections
- `exports/` — final assembled docs (if full draft)

## 5) Common warnings (what they mean)
- `Modal.Dict not available`: optional cache module not installed; safe to ignore.
- `429 signaled for semantic_scholar`: rate limit; use `SEMANTIC_SCHOLAR_API_KEY` or reduce parallelism.
- `Parallel query timeout`: one API timed out; reduce parallelism or raise timeouts.
- Gemini read timeout: set `ENABLE_GEMINI_GROUNDED=false` or increase `GEMINI_GROUNDED_TIMEOUT`.
- OpenAI-compatible read timeout: increase `OPENAI_TIMEOUT` and reduce `MODEL_MAX_OUTPUT_TOKENS`.

## 6) Notes on runtime
Deep research can be slow. Example observed run:
- Search-only (expose) with deep research cap 20 queries
- Runtime: ~14.5 minutes (13:32:00 → 13:46:33 on Jan 30, 2026)

Use this as a baseline for other runs.

