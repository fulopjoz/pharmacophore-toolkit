#!/usr/bin/env python3
"""AI-Researcher Integration Script for Pharmacophore-Toolkit Research.

This module provides utilities for integrating with the external AI-Researcher
tool for literature search, idea generation, and experiment planning.

The AI-Researcher tool (located at /home/dodo/Documents/projects/AI-Researcher/)
supports three modes:
1. Detailed Idea Description - Generate detailed plans from research ideas
2. Reference-Based Ideation - Generate new ideas from reference papers
3. Paper Generation Agent - Generate full paper drafts

Usage:
    # Start AI-Researcher first:
    # cd /home/dodo/Documents/projects/AI-Researcher && python web_ai_researcher.py
    
    # Then use this module:
    from experiments.ai_researcher_integration import AIResearcherClient
    
    client = AIResearcherClient()
    plan = client.generate_experiment_plan(
        idea="Implement grid-based spatial clustering for pharmacophore features",
        references=["PheSA 2024", "Consensus VS 2024"]
    )
"""

import json
import os
import requests
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Any


@dataclass
class ResearchQuery:
    """Research query template for AI-Researcher."""
    topic: str
    question: str
    context: str
    references: List[str]
    expected_output: str


# Predefined research queries for pharmacophore optimization
RESEARCH_QUERIES = {
    "grid_clustering": ResearchQuery(
        topic="Grid-Based Spatial Clustering for Molecules",
        question="What are the most efficient grid-based methods for clustering 3D molecular features?",
        context="""
        We have a pharmacophore toolkit that uses agglomerative hierarchical clustering
        to find consensus features from aligned molecules. Current complexity is O(n² log n).
        We want to implement faster grid-based clustering with O(n) complexity.
        
        Current implementation:
        - Feature extraction: [type, indices, x, y, z] format
        - Clustering: AgglomerativeClustering with distance_threshold
        - Parameters: tolerance (1.0-2.5 Å), occurrence_threshold (0.5-1.0)
        """,
        references=[
            "Clustering protein binding pockets (Stevenson 2023)",
            "PheSA pharmacophore-enhanced shape alignment (Wahl 2024)",
            "Spatial clustering with GNN (Pineda 2024)"
        ],
        expected_output="Implementation plan for grid-based clustering with benchmarks"
    ),
    
    "parameter_optimization": ResearchQuery(
        topic="Pharmacophore Parameter Optimization",
        question="How to automatically optimize tolerance and threshold for consensus pharmacophores?",
        context="""
        Our consensus pharmacophore models have two key parameters:
        - tolerance: spatial radius for feature clustering (default 2.0 Å)
        - occurrence_threshold: fraction of molecules that must share feature (default 0.5)
        
        Literature suggests:
        - Tolerance: 1.2-1.5 Å for flexibility
        - Threshold: 0.6-0.8 for partial matching
        
        We need automated parameter selection based on:
        - Target class (kinase, GPCR, etc.)
        - Number of reference molecules
        - Desired screening stringency
        """,
        references=[
            "Consensus holistic virtual screening (Moshawih 2024)",
            "Comparison of consensus pharmacophore screening (KEY_INSIGHTS.md)"
        ],
        expected_output="Active learning or Bayesian optimization approach for parameters"
    ),
    
    "gpu_acceleration": ResearchQuery(
        topic="GPU Acceleration for Molecular Distance Calculations",
        question="How to implement GPU-accelerated distance matrix computation for pharmacophore clustering?",
        context="""
        Current bottleneck: Computing pairwise distances between pharmacophore features.
        For N features, we compute N(N-1)/2 distances.
        
        With 100 molecules × 10 features = 1000 features → 500K distance calculations.
        
        Target: 10-100x speedup using CUDA or OpenCL.
        
        Constraints:
        - Must remain deterministic
        - Must work on standard consumer GPUs
        - Fallback to CPU for compatibility
        """,
        references=[
            "GPU acceleration in cheminformatics",
            "RAPIDS cuML for sklearn-compatible GPU clustering"
        ],
        expected_output="CUDA kernel design or cuML integration approach"
    ),
    
    "benchmark_datasets": ResearchQuery(
        topic="Benchmark Datasets for Pharmacophore Virtual Screening",
        question="Which benchmark datasets best evaluate pharmacophore-based virtual screening?",
        context="""
        We need to validate our consensus pharmacophore models on standard benchmarks.
        
        Known datasets:
        - DUD-E: 102 targets, ~20K actives, 50:1 decoy ratio
        - LIT-PCBA: 15 targets, higher quality actives
        - MUV: Maximum Unbiased Validation
        
        Metrics we compute: ROC-AUC, BEDROC, EF@1/5/10%, Youden's J
        
        Questions:
        - Which targets best test consensus methods?
        - How to handle chemotype diversity in actives?
        - What's acceptable performance (AUC > 0.7? 0.8?)
        """,
        references=[
            "DUD-E database paper",
            "LIT-PCBA benchmarking study",
            "Pharmacophore benchmark review papers"
        ],
        expected_output="Recommended benchmark protocol with target selection"
    )
}


class AIResearcherClient:
    """Client for interacting with AI-Researcher Gradio API."""
    
    def __init__(
        self,
        base_url: str = "http://127.0.0.1:7039",
        timeout: int = 300
    ):
        """Initialize client.
        
        Args:
            base_url: AI-Researcher Gradio server URL
            timeout: Request timeout in seconds
        """
        self.base_url = base_url
        self.timeout = timeout
    
    def health_check(self) -> bool:
        """Check if AI-Researcher server is running."""
        try:
            response = requests.get(f"{self.base_url}/", timeout=5)
            return response.status_code == 200
        except requests.RequestException:
            return False
    
    def generate_experiment_plan(
        self,
        idea: str,
        references: Optional[List[str]] = None,
        mode: str = "Detailed Idea Description"
    ) -> Dict[str, Any]:
        """Generate detailed experiment plan from research idea.
        
        Args:
            idea: Research idea description
            references: List of reference paper titles
            mode: AI-Researcher mode (Detailed Idea Description or Reference-Based Ideation)
        
        Returns:
            Dict with generated plan
            
        Note:
            This method requires the AI-Researcher Gradio server to be running.
            Start it with: cd AI-Researcher && python web_ai_researcher.py
        """
        if not self.health_check():
            return {
                "error": "AI-Researcher server not running",
                "suggestion": "Start with: cd AI-Researcher && python web_ai_researcher.py"
            }
        
        # Gradio API endpoint (standard for all Gradio apps)
        api_url = f"{self.base_url}/api/predict/"
        
        payload = {
            "data": [
                idea,
                "\n".join(references) if references else "",
                mode
            ]
        }
        
        try:
            response = requests.post(
                api_url,
                json=payload,
                timeout=self.timeout
            )
            
            if response.status_code == 200:
                return response.json()
            else:
                return {
                    "error": f"API error: {response.status_code}",
                    "response": response.text
                }
        except requests.RequestException as e:
            return {"error": str(e)}
    
    def run_research_query(self, query_name: str) -> Dict[str, Any]:
        """Run a predefined research query.
        
        Args:
            query_name: Key from RESEARCH_QUERIES dict
        
        Returns:
            API response dict
        """
        if query_name not in RESEARCH_QUERIES:
            available = list(RESEARCH_QUERIES.keys())
            return {"error": f"Unknown query. Available: {available}"}
        
        query = RESEARCH_QUERIES[query_name]
        
        idea = f"""
Topic: {query.topic}

Question: {query.question}

Context:
{query.context}

Expected Output: {query.expected_output}
"""
        
        return self.generate_experiment_plan(
            idea=idea,
            references=query.references,
            mode="Detailed Idea Description"
        )


def create_research_session(output_dir: Optional[Path] = None) -> Dict[str, Any]:
    """Create a research session with predefined queries.
    
    This function creates a session file that can be used to track
    research progress and results.
    
    Args:
        output_dir: Directory to save session file
    
    Returns:
        Session metadata dict
    """
    if output_dir is None:
        output_dir = Path(__file__).parent.parent / "docs" / "research" / "sessions"
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    from datetime import datetime
    session_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    session = {
        "session_id": session_id,
        "created": datetime.now().isoformat(),
        "queries": list(RESEARCH_QUERIES.keys()),
        "status": {q: "pending" for q in RESEARCH_QUERIES},
        "results": {}
    }
    
    session_file = output_dir / f"research_session_{session_id}.json"
    with open(session_file, 'w') as f:
        json.dump(session, f, indent=2)
    
    print(f"Research session created: {session_file}")
    return session


def main():
    """CLI for AI-Researcher integration."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="AI-Researcher integration for pharmacophore research"
    )
    parser.add_argument(
        '--check',
        action='store_true',
        help='Check if AI-Researcher server is running'
    )
    parser.add_argument(
        '--query',
        choices=list(RESEARCH_QUERIES.keys()),
        help='Run a predefined research query'
    )
    parser.add_argument(
        '--list-queries',
        action='store_true',
        help='List available research queries'
    )
    parser.add_argument(
        '--create-session',
        action='store_true',
        help='Create a new research session'
    )
    args = parser.parse_args()
    
    client = AIResearcherClient()
    
    if args.check:
        if client.health_check():
            print("✅ AI-Researcher server is running")
        else:
            print("❌ AI-Researcher server not running")
            print("   Start with: cd AI-Researcher && python web_ai_researcher.py")
        return
    
    if args.list_queries:
        print("\nAvailable Research Queries:")
        print("=" * 60)
        for name, query in RESEARCH_QUERIES.items():
            print(f"\n{name}:")
            print(f"  Topic: {query.topic}")
            print(f"  Question: {query.question[:60]}...")
        return
    
    if args.create_session:
        create_research_session()
        return
    
    if args.query:
        print(f"\nRunning query: {args.query}")
        print("=" * 60)
        
        query = RESEARCH_QUERIES[args.query]
        print(f"Topic: {query.topic}")
        print(f"Question: {query.question}")
        print("-" * 60)
        
        result = client.run_research_query(args.query)
        
        if "error" in result:
            print(f"\nError: {result['error']}")
            if "suggestion" in result:
                print(f"Suggestion: {result['suggestion']}")
        else:
            print("\nResult:")
            print(json.dumps(result, indent=2))
        return
    
    # Default: show help
    parser.print_help()


if __name__ == "__main__":
    main()
