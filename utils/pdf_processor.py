#!/usr/bin/env python3
"""Batch PDF text extraction utility using PyMuPDF.

Extracts text from PDFs page-by-page with configurable safety limits
to prevent memory/context overflow. Saves results as .extracted.md files.

Adapted from OpenDraft's pdf_fetcher.py dual-limit pattern.

Usage:
    python utils/pdf_processor.py --input-dir docs/research/papers/priority/
    python utils/pdf_processor.py --input-dir docs/ --max-pages 30 --recursive
"""

import argparse
import sys
from datetime import datetime
from pathlib import Path

try:
    import fitz  # PyMuPDF
except ImportError:
    print("ERROR: PyMuPDF not installed. Run: pip install PyMuPDF")
    sys.exit(1)


def extract_pdf_text(
    pdf_path: Path, max_pages: int = 20, max_chars: int = 50000
) -> dict:
    """Extract text from PDF with page + character limits.

    Returns dict with 'text', 'pages_extracted', 'total_pages', 'chars'.
    """
    try:
        doc = fitz.open(str(pdf_path))
        total_pages = doc.page_count
        text_parts = []
        pages_extracted = 0

        for i in range(min(max_pages, total_pages)):
            page = doc.load_page(i)
            text_parts.append(page.get_text("text"))
            pages_extracted += 1
            if sum(len(p) for p in text_parts) >= max_chars:
                break

        doc.close()
        text = "\n".join(text_parts)[:max_chars]
        return {
            "text": text,
            "pages_extracted": pages_extracted,
            "total_pages": total_pages,
            "chars": len(text),
        }
    except Exception as e:
        return {
            "text": "",
            "pages_extracted": 0,
            "total_pages": 0,
            "chars": 0,
            "error": str(e),
        }


def process_directory(
    input_dir: Path,
    output_dir: Path | None,
    max_pages: int,
    max_chars: int,
    recursive: bool,
    force: bool,
) -> list[dict]:
    """Process all PDFs in a directory. Returns list of result metadata."""
    if output_dir is None:
        output_dir = input_dir

    pattern = "**/*.pdf" if recursive else "*.pdf"
    pdfs = sorted(input_dir.glob(pattern))

    if not pdfs:
        print(f"No PDFs found in {input_dir}")
        return []

    print(f"Found {len(pdfs)} PDF(s) in {input_dir}")
    results = []

    for pdf_path in pdfs:
        rel = pdf_path.relative_to(input_dir)
        out_path = output_dir / rel.with_suffix(".extracted.md")
        out_path.parent.mkdir(parents=True, exist_ok=True)

        if out_path.exists() and not force:
            print(f"  SKIP (exists): {rel}")
            # Read existing metadata from first line
            first_line = out_path.read_text().split("\n")[0]
            results.append(
                {
                    "file": str(rel),
                    "output": str(out_path.relative_to(output_dir)),
                    "status": "skipped",
                    "size_mb": pdf_path.stat().st_size / (1024 * 1024),
                }
            )
            continue

        print(f"  Processing: {rel} ({pdf_path.stat().st_size / (1024*1024):.1f} MB)")
        result = extract_pdf_text(pdf_path, max_pages, max_chars)

        if result.get("error"):
            print(f"    ERROR: {result['error']}")
            results.append(
                {
                    "file": str(rel),
                    "status": "error",
                    "error": result["error"],
                    "size_mb": pdf_path.stat().st_size / (1024 * 1024),
                }
            )
            continue

        # Write extracted text as markdown
        header = (
            f"# {pdf_path.stem}\n\n"
            f"**Source**: `{rel}`\n"
            f"**Pages extracted**: {result['pages_extracted']}/{result['total_pages']}\n"
            f"**Characters**: {result['chars']:,}\n"
            f"**Extracted**: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n"
            f"---\n\n"
        )
        out_path.write_text(header + result["text"])
        print(
            f"    -> {out_path.name} "
            f"({result['pages_extracted']}/{result['total_pages']} pages, "
            f"{result['chars']:,} chars)"
        )

        results.append(
            {
                "file": str(rel),
                "output": str(out_path.relative_to(output_dir)),
                "status": "ok",
                "pages": f"{result['pages_extracted']}/{result['total_pages']}",
                "chars": result["chars"],
                "size_mb": pdf_path.stat().st_size / (1024 * 1024),
            }
        )

    # Generate index file
    index_path = output_dir / "_INDEX.md"
    write_index(index_path, results, input_dir)
    print(f"\nIndex written to: {index_path}")
    return results


def write_index(index_path: Path, results: list[dict], input_dir: Path):
    """Write summary index of all processed PDFs."""
    lines = [
        f"# PDF Extraction Index\n",
        f"**Directory**: `{input_dir}`\n",
        f"**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n",
        f"**Total files**: {len(results)}\n\n",
        "| File | Size | Pages | Chars | Status |\n",
        "|------|------|-------|-------|--------|\n",
    ]
    for r in results:
        size = f"{r.get('size_mb', 0):.1f} MB"
        pages = r.get("pages", "-")
        chars = f"{r.get('chars', 0):,}" if r.get("chars") else "-"
        status = r.get("status", "?")
        lines.append(f"| {r['file']} | {size} | {pages} | {chars} | {status} |\n")

    index_path.write_text("".join(lines))


def main():
    parser = argparse.ArgumentParser(
        description="Extract text from PDFs using PyMuPDF with safety limits"
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("docs/research/papers/priority"),
        help="Directory containing PDFs (default: docs/research/papers/priority)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory for .md files (default: same as input-dir)",
    )
    parser.add_argument(
        "--max-pages",
        type=int,
        default=20,
        help="Maximum pages to extract per PDF (default: 20)",
    )
    parser.add_argument(
        "--max-chars",
        type=int,
        default=50000,
        help="Maximum characters per PDF (default: 50000)",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Process subdirectories recursively",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing .extracted.md files",
    )

    args = parser.parse_args()

    if not args.input_dir.is_dir():
        print(f"ERROR: {args.input_dir} is not a directory")
        sys.exit(1)

    process_directory(
        args.input_dir,
        args.output_dir,
        args.max_pages,
        args.max_chars,
        args.recursive,
        args.force,
    )


if __name__ == "__main__":
    main()
