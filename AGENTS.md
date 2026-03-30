# Repository Guidelines

## Project Structure & Module Organization
This repository is currently a lightweight prototype. The top level contains [README.md](/workspaces/PCOA-prototype/README.md), which describes the goal: subsample ENA read data and generate PCoA plots quickly using `seqkit`, `deblur`, and `scikit-bio`.

As the project grows, keep code in clearly named top-level directories:
- `src/` for pipeline or analysis code
- `tests/` for automated checks
- `data/` for small sample inputs only; do not commit large raw datasets
- `results/` or `plots/` for generated outputs that are safe to version

## Build, Test, and Development Commands
There is no formal build system in the repository yet. Prefer simple, reproducible CLI workflows and document new commands in [README.md](/workspaces/PCOA-prototype/README.md) when you add them.

Useful baseline commands:
- `git status` checks local changes before editing or reviewing
- `python -m pytest` runs tests once a Python test suite exists
- `python script_name.py --help` verifies a script’s CLI contract

If you add a dependency manager such as `pip`, `uv`, or `conda`, include setup steps and one canonical local run command.

## Coding Style & Naming Conventions
Use 4-space indentation for Python. Prefer small, composable scripts or modules over notebook-only logic. Name Python files and functions in `snake_case`, classes in `PascalCase`, and keep command-line entry points descriptive, for example `subsample_reads.py`.

Write code that makes pipeline stages obvious: input loading, subsampling, QC, distance calculation, and plotting should be easy to trace.

## Testing Guidelines
Add tests under `tests/` with names like `test_subsampling.py`. Focus first on deterministic units such as argument parsing, sample filtering, and distance-matrix generation. For data-heavy steps, use tiny fixtures rather than full ENA inputs.

Run tests locally before opening a PR. If automated coverage is introduced, document the expected threshold in this file and the README.

## Commit & Pull Request Guidelines

- Delete unused or obsolete files when your changes make them irrelevant (refactors, feature removals, etc.), and revert files only when the change is yours or explicitly requested. If a git operation leaves you unsure about other agents' in-flight work, stop and coordinate instead of deleting.
- **Before attempting to delete a file to resolve a local type/lint failure, stop and ask the user.** Other agents are often editing adjacent files; deleting their work to silence an error is never acceptable without explicit approval.
- NEVER edit `.env` or any environment variable files—only the user may change them.
- Coordinate with other agents before removing their in-progress edits—don't revert or delete work you didn't author unless everyone agrees.
- Moving/renaming and restoring files is allowed.
- ABSOLUTELY NEVER run destructive git operations (e.g., `git reset --hard`, `rm`, `git checkout`/`git restore` to an older commit) unless the user gives an explicit, written instruction in this conversation. Treat these commands as catastrophic; if you are even slightly unsure, stop and ask before touching them. *(When working within Cursor or Codex Web, these git limitations do not apply; use the tooling's capabilities as needed.)*
- Never use `git restore` (or similar commands) to revert files you didn't author—coordinate with other agents instead so their in-progress work stays intact.
- Always double-check git status before any commit
- Keep commits atomic: commit only the files you touched and list each path explicitly. For tracked files run `git commit -m "<scoped message>" -- path/to/file1 path/to/file2`. For brand-new files, use the one-liner `git restore --staged :/ && git add "path/to/file1" "path/to/file2" && git commit -m "<scoped message>" -- path/to/file1 path/to/file2`.
- Quote any git paths containing brackets or parentheses (e.g., `src/app/[candidate]/**`) when staging or committing so the shell does not treat them as globs or subshells.
- When running `git rebase`, avoid opening editors—export `GIT_EDITOR=:` and `GIT_SEQUENCE_EDITOR=:` (or pass `--no-edit`) so the default messages are used automatically.
- Never amend commits unless you have explicit written approval in the task thread.
