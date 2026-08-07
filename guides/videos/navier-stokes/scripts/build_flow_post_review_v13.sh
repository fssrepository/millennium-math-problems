#!/usr/bin/env bash
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
guide_dir="$(cd -- "${script_dir}/.." && pwd)"
video_input="${1:-/home/raxim/Downloads/Initial_Scene_-_2026-08-05_202608051411.mp4}"
audio_input="${2:-${guide_dir}/renders/edit/navier-stokes-flow-post-v09-review-720p.mp4}"
output_file="${3:-${guide_dir}/renders/edit/navier-stokes-flow-post-v13-review-720p.mp4}"
exec "${script_dir}/build_flow_post_review_v12.sh" "${video_input}" "${audio_input}" "${output_file}"
