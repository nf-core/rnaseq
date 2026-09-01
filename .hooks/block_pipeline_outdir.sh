#!/usr/bin/env bash
# This hook is used to block commits if they include staged files inside a directory
# which also contains a subdirectory called `pipeline_info`. The purpose of this is to
# prevent users from inadvertently committing output from pipeline test runs inside the
# development directory.

set -e

status=0
seen_dirs=""

while IFS= read -r file; do
    # The offending output bundle's root is the ancestor directory that has
    # `pipeline_info` as an immediate child, so callers can restore it in one go.
    if [[ "$file" == pipeline_info/* ]]; then
        top_dir="pipeline_info"
    elif [[ "$file" == */pipeline_info/* ]]; then
        top_dir="${file%%/pipeline_info/*}"
    else
        top_dir=""
        dir=$(dirname "$file")
        while [[ "$dir" != "." && "$dir" != "/" ]]; do
            if [[ -d "$dir/pipeline_info" ]]; then
                top_dir="$dir"
                break
            fi
            dir=$(dirname "$dir")
        done
    fi

    if [[ -n "$top_dir" ]]; then
        echo "❌ Commit blocked: Please do not commit output from pipeline test runs to the pipeline code itself: $file"
        status=1
        case "$seen_dirs" in
            *"|$top_dir|"*) ;;
            *)
                echo "Run 'git restore --staged $top_dir' to remove the whole output folder from the staging area."
                seen_dirs="$seen_dirs|$top_dir|"
                ;;
        esac
    fi
done < <(git diff --cached --name-only)

exit "$status"
