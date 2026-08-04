#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
guide_dir="$(cd -- "${script_dir}/.." && pwd)"
input_file="${1:-${guide_dir}/renders/flow/drafts/S05_flow_lite_720p_take01.mp4}"
output_file="${2:-${guide_dir}/renders/flow/edits/S05_flow_lite_take01_3p83s_clean_silent.mp4}"

if [[ ! -f "${input_file}" ]]; then
  echo "Missing input: ${input_file}" >&2
  exit 1
fi

mkdir -p -- "$(dirname -- "${output_file}")"

# The generated scan mutates into a triangular beam with a circular handle from
# roughly 1.8 to 5.8 seconds. Use the clean 6.5--8.0 second fluid tail, slow it
# to the final S05 duration. The malformed generated scan is not replaced in
# the base plate: the final evidence marker and pointer are added with the
# exact typography pass. Source audio is intentionally discarded because
# slowing it created audible artifacts; final sound is assembled separately.
ffmpeg -loglevel error -y \
  -i "${input_file}" \
  -filter_complex "[0:v]trim=start=6.5:end=8,setpts=(PTS-STARTPTS)*2.553333333,fps=30,trim=duration=3.83,setpts=PTS-STARTPTS,format=yuv420p[outv]" \
  -map "[outv]" \
  -t 3.83 \
  -c:v libx264 \
  -preset medium \
  -crf 18 \
  -an \
  -movflags +faststart \
  "${output_file}"

ffprobe -v error \
  -show_entries format=duration,size:stream=index,codec_type,codec_name,width,height,pix_fmt,r_frame_rate \
  -of default=noprint_wrappers=1 \
  "${output_file}"
