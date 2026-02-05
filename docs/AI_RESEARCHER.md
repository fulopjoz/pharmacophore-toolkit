# AI-Researcher (HKUDS) Setup + e-infra Models

This repo is set up to use a local clone of `HKUDS/AI-Researcher` located at `/home/dodo/Documents/projects/AI-Researcher`.

## 1) One-time install

```bash
cd /home/dodo/Documents/projects/AI-Researcher
python -m venv .venv
. .venv/bin/activate
python -m pip install -U pip
python -m pip install -e .
playwright install chromium
```

## 2) Configure e-infra (OpenAI-compatible)

Copy the provided example and fill your key:

```bash
cd /home/dodo/Documents/projects/AI-Researcher
cp .env.einfra.example .env
```

Recommended model mapping for `https://llm.ai.e-infra.cz/v1`:

- `COMPLETION_MODEL=openai/deepseek-v3.2-thinking` (best reasoning/planning)
- `CHEEP_MODEL=openai/deepseek-v3.2` (cheaper helper model)
- `EMBEDDING_MODEL=qwen3-embedding-4b` (RAG embeddings)
- `LITELLM_PROVIDER=openai` (ensures LiteLLM provider is set for custom base URL)
- `BASE_IMAGES=tjbtech1/airesearcher:v1` (Docker image used by the agent)
- `GPUS=` (leave blank to run CPU-only; set `"device=0"` if you want GPU)
- `GITHUB_AI_TOKEN=` (optional; without it GitHub search is skipped)

## 3) Quick connectivity test (no repo-specific logic)

```bash
cd /home/dodo/Documents/projects/AI-Researcher
. .venv/bin/activate
python smoke_test_llm.py
```
`smoke_test_llm.py` auto-loads `.env` if present.

## 4) Run the Web UI

```bash
cd /home/dodo/Documents/projects/AI-Researcher
. .venv/bin/activate
python web_ai_researcher.py
```
If port 7039 is in use, set `GRADIO_SERVER_PORT=7040` (or any free port) before launch.
If you see `tjbtech1/paperapp:latest` errors, set `BASE_IMAGES=tjbtech1/airesearcher:v1` in `.env`.

## 5) Using AI-Researcher from another Python project

Practical integration (recommended): run it as a separate process from your project so it can manage its own working directory and `.env`.

Example (subprocess):

```python
import os
import subprocess

AI_RESEARCHER_DIR = "/home/dodo/Documents/projects/AI-Researcher"

env = os.environ.copy()
env["API_BASE_URL"] = "https://llm.ai.e-infra.cz/v1"
env["OPENAI_API_KEY"] = os.environ["OPENAI_API_KEY"]
env["COMPLETION_MODEL"] = "deepseek-v3.2-thinking"
env["CHEEP_MODEL"] = "deepseek-v3.2"
env["EMBEDDING_MODEL"] = "qwen3-embedding-4b"

subprocess.run(
    [f"{AI_RESEARCHER_DIR}/.venv/bin/python", f"{AI_RESEARCHER_DIR}/web_ai_researcher.py"],
    env=env,
    check=True,
)
```

Notes:
- Full “implement + train in container” flows require Docker + the upstream image (`tjbtech1/airesearcher:v1`). If Docker isn’t installed, only the non-container parts will work.
