#!/bin/bash
#
# Apply one patch spanning several contribs and prepare patch releases.
#
# Usage:
#   scripts/internal/apply-multi-contrib-patch.sh <patch-file>
#
# Run this script from the top-level fjcontrib directory. Each contrib
# affected by the patch must initially be at the release listed in
# contribs.svn, and its trunk must have no unreleased content changes.

# Load common setup and utilities.
script_dir=$(cd "$(dirname "$0")" && pwd)
. "$script_dir/common.sh"

top_dir=$(cd "$script_dir/../.." && pwd)
if [[ "$PWD" != "$top_dir" ]]; then
  echo "Run this script from the top-level fjcontrib directory:"
  echo "  cd $top_dir"
  echo "  scripts/internal/$(basename "$0") <patch-file>"
  exit 1
fi

if [[ $# -eq 1 && ("$1" == "-h" || "$1" == "--help") ]]; then
  echo "Usage:"
  echo "  scripts/internal/$(basename "$0") <patch-file>"
  exit 0
fi

if [[ $# -ne 1 ]]; then
  echo "Usage:"
  echo "  scripts/internal/$(basename "$0") <patch-file>"
  exit 1
fi

if [[ ! -f "$1" || ! -r "$1" ]]; then
  echo "Cannot read patch file: $1"
  exit 1
fi

patch_dir=$(cd "$(dirname "$1")" && pwd)
patch_file="$patch_dir/$(basename "$1")"
if [[ ! -s "$patch_file" ]]; then
  echo "Patch file is empty: $patch_file"
  exit 1
fi

message_file=""
diagnostics_file=""
cleanup() {
  if [[ -n "$message_file" && -e "$message_file" ]]; then
    rm -f "$message_file"
  fi
  if [[ -n "$diagnostics_file" && -e "$diagnostics_file" ]]; then
    rm -f "$diagnostics_file"
  fi
}
trap cleanup EXIT
trap 'cleanup; exit 1' HUP INT TERM

# Extract all paths mentioned in common unified, SVN and Git patch formats.
# Contrib filenames do not contain whitespace, so paths can safely be kept in
# whitespace-separated lists, as is done elsewhere in the fjcontrib scripts.
patch_paths=$(awk '
  $1 == "diff" && $2 == "--git" {
    print $3
    print $4
    next
  }
  $1 == "Index:" {
    print $2
    next
  }
  $1 == "---" || $1 == "+++" {
    print $2
  }
' "$patch_file" | sed '/^\/dev\/null$/d;s|^[ab]/||' | sort -u)

if [[ -z "$patch_paths" ]]; then
  echo "No file paths could be identified in $patch_file"
  exit 1
fi

contribs=""
for path in $patch_paths; do
  if [[ "$path" == /* || "$path" == ".." || "$path" == ../* ||
        "$path" == */../* || "$path" != */* ]]; then
    echo "Unsafe or non-contrib path in patch: $path"
    exit 1
  fi

  contrib=${path%%/*}
  get_contrib_version "$contrib" contribs.svn release_location
  if [[ "$release_location" == "[None]" ]]; then
    echo "Patch path does not belong to a contrib in contribs.svn: $path"
    exit 1
  fi
  if ! item_is_in_list "$contrib" "$contribs"; then
    contribs="$contribs $contrib"
  fi
done
contribs=${contribs# }

if [[ -z "$contribs" ]]; then
  echo "The patch does not affect any contribs"
  exit 1
fi

read_local_version() {
  local contrib_dir=$1
  local value

  if [[ -f "$contrib_dir/FJCONTRIB.cfg" ]]; then
    value=$(sed -n 's/^[[:space:]]*version[[:space:]]*:[[:space:]]*//p' \
      "$contrib_dir/FJCONTRIB.cfg" | head -n1)
  elif [[ -f "$contrib_dir/VERSION" ]]; then
    value=$(head -n1 "$contrib_dir/VERSION")
  else
    return 1
  fi

  [[ -n "$value" ]] || return 1
  printf '%s\n' "$value"
}

increment_patch_version() {
  local version=$1
  local major
  local minor
  local patch

  major=${version%%.*}
  version=${version#*.}
  minor=${version%%.*}
  patch=${version#*.}
  printf '%s.%s.%s\n' "$major" "$minor" "$((patch + 1))"
}

check_trunk_matches_release() {
  local contrib=$1
  local release_location=$2
  local release_url="$svn_read/contribs/$contrib/$release_location"
  local trunk_url="$svn_read/contribs/$contrib/trunk"
  local raw_diff
  local meaningful_diff

  raw_diff=$(mktemp "${TMPDIR:-/tmp}/fjcontrib-trunk-diff.XXXXXX") || return 1
  meaningful_diff=$(mktemp \
    "${TMPDIR:-/tmp}/fjcontrib-trunk-meaningful.XXXXXX") || {
    rm -f "$raw_diff"
    return 1
  }

  # -B ignores added or removed blank lines, while --strip-trailing-cr
  # ignores CRLF-versus-LF changes. SVN properties are ignored separately.
  if ! svn diff --ignore-properties --diff-cmd diff \
      -x '-u -B --strip-trailing-cr' \
      "$release_url" "$trunk_url" > "$raw_diff"; then
    echo "Could not compare $release_location and trunk for $contrib"
    rm -f "$raw_diff" "$meaningful_diff"
    return 1
  fi

  # SVN emits an Index heading even when the external diff suppresses all
  # changes. Remove that scaffolding before testing whether content differs.
  if ! sed '/^Index: /d;/^=\{3,\}$/d;/^[[:space:]]*$/d' \
      "$raw_diff" > "$meaningful_diff"; then
    echo "Could not process the release-to-trunk diff for $contrib"
    rm -f "$raw_diff" "$meaningful_diff"
    return 1
  fi

  if [[ -s "$meaningful_diff" ]]; then
    echo "Unreleased trunk changes found for $contrib:"
    sed 's/^/  /' "$meaningful_diff"
    rm -f "$raw_diff" "$meaningful_diff"
    return 1
  fi

  rm -f "$raw_diff" "$meaningful_diff"
  return 0
}

echo
echo "Checking affected contribs:"
echo "  $contribs"

index=0
for contrib in $contribs; do
  get_contrib_version "$contrib" contribs.svn release_location
  if [[ "$release_location" != tags/* ]]; then
    echo "$contrib is not at a released version in contribs.svn: $release_location"
    exit 1
  fi
  release_version=${release_location#tags/}
  if ! printf '%s\n' "$release_version" |
      grep -Eq '^(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)$'; then
    echo "$contrib has a non-semantic release version: $release_version"
    exit 1
  fi

  get_svn_info "$contrib" local_mode local_location
  if [[ "$local_location" != "$release_location" ]]; then
    echo "$contrib must initially be at $release_location"
    echo "  Current working copy: $local_location"
    exit 1
  fi
  if ! local_version=$(read_local_version "$contrib"); then
    echo "Could not read the version declared by the local $contrib release"
    exit 1
  fi
  if [[ "$local_version" != "$release_version" ]]; then
    echo "Version mismatch in the local release of $contrib:"
    echo "  contribs.svn:     $release_version"
    echo "  declared version: $local_version"
    exit 1
  fi
  check_pending_modifications "$contrib" || {
    echo "$contrib has pending modifications; clean or commit them first."
    exit 1
  }

  echo "  Comparing $release_location with trunk for $contrib"
  if ! check_trunk_matches_release "$contrib" "$release_location"; then
    exit 1
  fi

  old_versions[$index]=$release_version
  new_versions[$index]=$(increment_patch_version "$release_version")
  contrib_array[$index]=$contrib
  index=$((index + 1))
done

echo
echo "Switching affected contribs to trunk:"
for contrib in $contribs; do
  echo "$contrib:"
  "$script_dir/switch-to-version.sh" "$contrib" trunk || {
    echo "Failed to switch $contrib to trunk."
    exit 1
  }

  get_svn_info "$contrib" local_mode local_location
  if [[ "$local_location" != "trunk" ]]; then
    echo "$contrib did not switch to trunk (found $local_location)"
    exit 1
  fi
  if [[ ! -f "$contrib/NEWS" ]]; then
    echo "$contrib has no NEWS file"
    exit 1
  fi
done

# Git-format patches normally need one leading path component removed.
strip_level=0
if grep -Eq '^diff --git a/[^[:space:]]+ b/[^[:space:]]+' "$patch_file" ||
    { grep -Eq '^--- a/[^[:space:]]+' "$patch_file" &&
      grep -Eq '^\+\+\+ b/[^[:space:]]+' "$patch_file"; }; then
  strip_level=1
fi

diagnostics_file=$(mktemp "${TMPDIR:-/tmp}/fjcontrib-patch.XXXXXX") || exit 1
if ! patch -p"$strip_level" --dry-run < "$patch_file" > "$diagnostics_file" 2>&1; then
  alternate_strip_level=$((1 - strip_level))
  if patch -p"$alternate_strip_level" --dry-run < "$patch_file" \
      > "$diagnostics_file" 2>&1; then
    strip_level=$alternate_strip_level
  else
    echo
    echo "The patch does not apply cleanly to the contrib trunks:"
    sed 's/^/  /' "$diagnostics_file"
    exit 1
  fi
fi

message_file=$(mktemp "${TMPDIR:-/tmp}/fjcontrib-message.XXXXXX") || exit 1
{
  echo "# Replace XYZ and edit this into the one-sentence NEWS/ChangeLog entry."
  echo "# Lines beginning with # are ignored. Keep the entry on one line."
  echo "Batch application of patch from XYZ to resolve compilation issues with C++23"
} > "$message_file"

editor_command=${VISUAL:-${EDITOR:-vi}}
if ! $editor_command "$message_file"; then
  echo "Editor failed: $editor_command"
  exit 1
fi

message=$(sed '/^[[:space:]]*#/d;/^[[:space:]]*$/d' "$message_file")
message_line_count=$(printf '%s\n' "$message" | awk 'END { print NR }')
if [[ -z "$message" || "$message_line_count" -ne 1 ]]; then
  echo "The release note must be exactly one non-comment, non-empty line."
  exit 1
fi
if printf '%s\n' "$message" | grep -q 'XYZ'; then
  echo "The release note still contains the XYZ placeholder."
  exit 1
fi

echo
echo "Summary"
echo "-------"
echo "Patch:       $patch_file"
echo "Strip level: -p$strip_level"
echo "Release note: $message"
echo
printf "  %-32s %-12s %-12s\n" "Contrib" "Current" "New"
printf "  %-32s %-12s %-12s\n" "-------" "-------" "---"
for ((index = 0; index < ${#contrib_array[@]}; index++)); do
  printf "  %-32s %-12s %-12s\n" \
    "${contrib_array[$index]}" \
    "${old_versions[$index]}" \
    "${new_versions[$index]}"
done
echo
echo "Each trunk matches its release tag, ignoring blank lines, line endings,"
echo "and SVN properties. All contribs are now on trunk, and the patch passes"
echo "a dry run."
get_yesno_answer "Apply the patch and prepare these patch releases?"
if [[ $? -eq 0 ]]; then
  echo "Aborting without applying the patch or changing release metadata."
  echo "The affected contribs remain switched to trunk."
  exit 1
fi

echo
echo "Applying $patch_file"
patch -p"$strip_level" < "$patch_file" || {
  echo "Patch application failed. Inspect the working copies before continuing."
  exit 1
}

# Schedule files created or removed by the patch so the printed commit command
# includes them. Only paths explicitly named by the patch are considered.
for path in $patch_paths; do
  status_line=$(svn status "$path" 2>/dev/null | head -n1)
  status_code=${status_line:0:1}
  if [[ "$status_code" == "?" ]]; then
    svn add "$path" || exit 1
  elif [[ "$status_code" == "!" ]]; then
    svn rm --force "$path" || exit 1
  fi
done

update_version_file() {
  local contrib=$1
  local new_version=$2
  local version_file
  local temporary_file
  local updated_version

  if [[ -f "$contrib/FJCONTRIB.cfg" ]]; then
    version_file="$contrib/FJCONTRIB.cfg"
    temporary_file=$(mktemp "${TMPDIR:-/tmp}/fjcontrib-version.XXXXXX") || return 1
    sed -E \
      "s/^([[:space:]]*version[[:space:]]*:[[:space:]]*).*/\\1$new_version/" \
      "$version_file" > "$temporary_file" || return 1
  elif [[ -f "$contrib/VERSION" ]]; then
    version_file="$contrib/VERSION"
    temporary_file=$(mktemp "${TMPDIR:-/tmp}/fjcontrib-version.XXXXXX") || return 1
    awk -v version="$new_version" 'NR == 1 { print version; next } { print }' \
      "$version_file" > "$temporary_file" || return 1
  else
    echo "$contrib has neither FJCONTRIB.cfg nor VERSION"
    return 1
  fi

  cp "$temporary_file" "$version_file" || return 1
  rm -f "$temporary_file"

  if ! updated_version=$(read_local_version "$contrib"); then
    return 1
  fi
  [[ "$updated_version" == "$new_version" ]]
}

prepend_release_entries() {
  local contrib=$1
  local new_version=$2
  local metadata_file
  local temporary_file
  local path
  local relative_path

  if [[ -f "$contrib/FJCONTRIB.cfg" ]]; then
    metadata_file="FJCONTRIB.cfg (version)"
  else
    metadata_file="VERSION"
  fi

  temporary_file=$(mktemp "${TMPDIR:-/tmp}/fjcontrib-news.XXXXXX") || return 1
  {
    printf '%s: release of version %s\n' "$(date '+%Y/%m/%d')" "$new_version"
    printf '%s\n' "$message" | fold -s -w 76 | sed 's/^/            /'
    echo
    cat "$contrib/NEWS"
  } > "$temporary_file" || return 1
  cp "$temporary_file" "$contrib/NEWS" || return 1
  rm -f "$temporary_file"

  temporary_file=$(mktemp "${TMPDIR:-/tmp}/fjcontrib-changelog.XXXXXX") || return 1
  {
    printf '%s  %s\n\n' "$(date '+%Y-%m-%d')" "$changelog_author"
    printf '\t* NEWS:\n'
    printf '\t* %s:\n' "$metadata_file"
    printf '\trelease of version %s\n' "$new_version"
    echo
    for path in $patch_paths; do
      if [[ "$path" == "$contrib/"* ]]; then
        relative_path=${path#*/}
        printf '\t* %s:\n' "$relative_path"
      fi
    done
    printf '%s\n' "$message" | fold -s -w 72 | sed 's/^/\t/'
    echo
    if [[ -f "$contrib/ChangeLog" ]]; then
      cat "$contrib/ChangeLog"
    fi
  } > "$temporary_file" || return 1
  cp "$temporary_file" "$contrib/ChangeLog" || return 1
  rm -f "$temporary_file"

  status_line=$(svn status "$contrib/ChangeLog" 2>/dev/null | head -n1)
  if [[ "${status_line:0:1}" == "?" ]]; then
    svn add "$contrib/ChangeLog" || return 1
  fi
}

changelog_author=${CHANGELOG_AUTHOR:-$(id -F 2>/dev/null)}
changelog_author=${changelog_author:-$(whoami)}

for ((index = 0; index < ${#contrib_array[@]}; index++)); do
  contrib=${contrib_array[$index]}
  new_version=${new_versions[$index]}
  update_version_file "$contrib" "$new_version" || {
    echo "Failed to update the version for $contrib"
    exit 1
  }
  prepend_release_entries "$contrib" "$new_version" || {
    echo "Failed to update NEWS and ChangeLog for $contrib"
    exit 1
  }
done

echo
echo "Patch and release metadata updates are complete."
echo "Review the diffs and run the appropriate tests. Then commit each contrib:"
echo
for ((index = 0; index < ${#contrib_array[@]}; index++)); do
  printf '  cd %q; svn commit -m %q; cd ..\n' \
    "${contrib_array[$index]}" "$message"
done
echo
echo "After all commits succeed, release each contrib:"
echo
for contrib in "${contrib_array[@]}"; do
  echo "  scripts/release-contrib.sh $contrib"
done

echo
get_yesno_answer "Run the commit and release commands above now?"
if [[ $? -eq 1 ]]; then
  echo
  echo "Committing each contrib:"
  for ((index = 0; index < ${#contrib_array[@]}; index++)); do
    contrib=${contrib_array[$index]}
    echo "  svn commit $contrib"
    svn commit -m "$message" "$contrib" || {
      echo "Failed to commit $contrib. Remaining commits and releases were not run."
      exit 1
    }
  done

  echo
  echo "Releasing each contrib:"
  for contrib in "${contrib_array[@]}"; do
    echo "  scripts/release-contrib.sh $contrib"
    "$script_dir/../release-contrib.sh" "$contrib" || {
      echo "Failed to release $contrib. Remaining releases were not run."
      exit 1
    }
  done
else
  echo "Leaving the prepared changes for manual review and execution."
fi
