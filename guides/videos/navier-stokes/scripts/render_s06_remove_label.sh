#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
guide_dir="$(cd -- "${script_dir}/.." && pwd)"
input_file="${1:-${guide_dir}/renders/flow/drafts/S06_flow_lite_720p_take01_rejected_text_and_markers.mp4}"
output_file="${2:-${guide_dir}/renders/flow/edits/S06_flow_lite_take01_8s_label_removed.mp4}"

if [[ ! -f "${input_file}" ]]; then
  echo "Missing input: ${input_file}" >&2
  exit 1
fi

mkdir -p -- "$(dirname -- "${output_file}")"

# Copy a small dark background patch over the generated 4%/6% label. Preserve
# the complete eight-second timeline, all generated motion and source audio.
ffmpeg -loglevel error -y \
  -i "${input_file}" \
  -filter_complex "[0:v]split=2[base][patch];[patch]crop=136:88:90:320,boxblur=4:2[clean];[base][clean]overlay=240:320[outv]" \
  -map "[outv]" -map 0:a? \
  -c:v libx264 -preset medium -crf 18 -pix_fmt yuv420p \
  -c:a copy -movflags +faststart \
  "${output_file}"

ffprobe -v error \
  -show_entries format=duration,size:stream=index,codec_type,codec_name,width,height,r_frame_rate,sample_rate,channels \
  -of default=noprint_wrappers=1 \
  "${output_file}"
