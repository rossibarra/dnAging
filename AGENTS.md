1. Git repo check
- If the current directory is not a git repository, ask the user whether they want to create one before making changes.

2. Default permission when the file is recoverable
- Permission to modify an existing data or markdown file is granted by default when the file is made **recoverable** first. A file is recoverable when either:
  - it is tracked by git and **clean since the last commit** (so `git checkout -- <path>` restores it), or
  - you preserve its current contents in a `.bak` file **before** modifying it, following rule 6.
- If the file is untracked or already has uncommitted changes, create or append to its `.bak` automatically and proceed without asking.
- Permission is not required for creating a new data file.

3. Data-file definition
- A data file is any non-code project file, especially plain-text files used as inputs, outputs, configuration, or reference data.
- Examples include: `*.txt`, `*.csv`, `*.tsv`, `*.fastq`, `*.fastq.gz`, `*.sam`, `*.maf`, `*.yaml`, `*.yml`, `*.json`, `*.map`, `*.bed`, `*.vcf`, `*.gvcf`, `*.fai`, and similar tabular or reference files.
- **Markdown (`*.md`) is documentation, not data** (README, MATH, NOTES and similar prose files).
- The automatic permission and recoverability requirements in rule 2 apply equally to markdown.
- If unsure whether a file is a data file, treat it as a data file and ask permission first.

4. Uncommitted changes check (target file only)
- Before modifying a file, check whether the **target file** has uncommitted changes.
- If the target file has uncommitted changes, preserve its current contents under rule 6 and proceed without asking. Do not commit the user's changes unless the user explicitly requests a commit.
- This rule is target-file-only (not repo-wide).

5. One permission can cover multiple files
- A single permission request is sufficient if the user clearly authorizes modifying multiple specific data files in the same task.

6. Backups for recoverability
- A `.bak` snapshot is required when an existing target file is **not** tracked-and-clean. Preserve it before modifying, using the `.bak` suffix (for example, `file.txt.bak`).
- If no `.bak` exists, copy the current file to `<file>.bak`.
- If `<file>.bak` already exists, do not overwrite it. Append a timestamped snapshot of the target file's current contents to the existing backup, with a clear delimiter containing an ISO-8601 timestamp and the source path.
- When the file is tracked by git and clean since the last commit, git already holds the pristine copy and no `.bak` is needed.
- `*.bak` is gitignored, so backups never enter a commit.

7. Symlink write policy
- If a path to be modified is a symbolic link, never modify the symlink target.
- If modification is required (with permission), create a regular-file copy in the current working directory, modify that copy, and replace the symlink path with the modified regular file.
- Prefer `path.tmp` + atomic rename (`mv path.tmp path`) so the symlink is replaced by a regular file.
- Do not use in-place editors (`sed -i`, `perl -pi`, etc.) on symlink paths.
- If a backup is required, back up the symlink path as it exists before replacement.
- If a symlink is replaced, create or append an entry in `symlinks.md` (in the repository root) recording: the original symlink path, the original symlink target path, and the replacement file path (the path after replacement).

8. Write scope restriction
- Do not create, modify, or delete files outside the current working directory unless the user explicitly requests it.
- Exception: temporary files may be created or modified in `/tmp` and system temporary directories (for example, macOS `/var/folders/...`) when needed for task execution.
- Files written to `/tmp` should be treated as temporary working files, not final outputs, unless the user explicitly requests otherwise.
- This rule does not permit modifying symlink targets outside the current working directory; symlink paths must follow the symlink write policy.
- Tool-generated temporary files, caches, and logs are allowed in the current working directory or approved temporary directories when required to complete the task.


9. Data-file permission exception (temporary paths)
- Permission is not required for creating or modifying data files under `/tmp` or system temporary directories (for example macOS `/var/folders/...`) when they are temporary working files used to run, test, or validate the project.
- This exception does not apply to files in the repository working directory (including `tests/`, `results/`, `example_data/`, `config.yaml`, `README.md`, `*.md`, etc.) or any other non-temporary location.
- Temporary files created under `/tmp` remain subject to the symlink write policy if the target path is a symlink.
- The agent should prefer `/tmp/<project>-<purpose>/...` paths for temporary data outputs to make scope explicit.

10. Bootstrap AGENTS.md into repos

If AGENTS.md is missing in the current git repository root, copy the default AGENTS.md into that repository before doing other work.
Do not overwrite an existing AGENTS.md.
If the current directory is not a git repository, follow the existing git-repo check rule first.

11. Conda environment bootstrap

- Before running any Python, pytest, pip, or project CLI command, first initialize conda in the shell.
- If the `module` command is available, use `module load conda`.
- If `module` is not available but `conda` is installed, initialize conda by sourcing `$(conda info --base)/etc/profile.d/conda.sh`.
- If the current working directory contains `environment.yml`, detect the environment name from its `name:` field and run `conda activate <that-name>` before continuing.
- If both `environment.yml` and an already-active conda environment are present, prefer the environment named in `environment.yml`.
- Run Python-related commands in a login bash shell so `module` and `conda activate` work correctly.
- If neither `module load conda` nor direct `conda` initialization is available, stop and report the error before running project commands in another Python environment.
- If activation fails, stop and report the error before running project commands in another Python environment.

12. No routine permission prompt for recoverable edits
- Do not ask for routine permission to modify an existing data or markdown file when rule 2 can make it recoverable automatically.
- A commit is not part of the automatic backup workflow. Commit only when the user explicitly requests it.
- Ask only when safe recoverability cannot be established or another rule independently requires confirmation.
