#!/bin/bash
# ==============================================================================
# setup_github.sh
# Run this once to initialize the repository and push to GitHub
# ==============================================================================

# 1. Initialize git
git init
git add README.md ROADMAP.md .gitignore R/00_functions.R

# Create empty placeholder files so empty directories are tracked
touch data/.gitkeep
touch output/tables/.gitkeep
touch output/figures/.gitkeep
touch output/notes/.gitkeep
touch tex/.gitkeep

git add data/.gitkeep output/tables/.gitkeep output/figures/.gitkeep \
        output/notes/.gitkeep tex/.gitkeep

# 2. First commit
git commit -m "Initial commit: repo structure, README, ROADMAP, core functions

- README.md: project overview, data instructions, run guide
- ROADMAP.md: reading list (SUW 2025 done), thesis outline, task tracker
- R/00_functions.R: modular kappa weight functions (tau_u, tau_a10, unnorm)
- .gitignore: excludes data files and R artifacts"

# 3. Connect to GitHub (run after creating repo on github.com)
# Replace with your actual GitHub username and repo name
echo ""
echo "Now go to https://github.com/new and create a repo called:"
echo "  thesis-kappa-late"
echo ""
echo "Then run:"
echo "  git remote add origin https://github.com/YOUR_USERNAME/thesis-kappa-late.git"
echo "  git branch -M main"
echo "  git push -u origin main"
