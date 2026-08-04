#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
guide_dir="$(cd -- "${script_dir}/.." && pwd)"
input_file="${1:-${guide_dir}/renders/flow/drafts/S02_flow_lite_720p_take02.mp4}"
overlay_file="${2:-${guide_dir}/assets/overlay-top-gradient-1080x1920-rgba.png}"
output_file="${3:-${guide_dir}/renders/flow/edits/S02_flow_lite_take02_6s_gradient_freeze.mp4}"

for required_file in "${input_file}" "${overlay_file}"; do
  if [[ ! -f "${required_file}" ]]; then
    echo "Missing input: ${required_file}" >&2
    exit 1
  fi
done

mkdir -p -- "$(dirname -- "${output_file}")"

# Keep approved motion through 5.375 seconds, bake the subtle upper-screen
# text-safe gradient, then hold the last approved frame for 0.625 seconds.
ffmpeg -loglevel error -y \
  -i "${input_file}" \
  -loop 1 -i "${overlay_file}" \
  -filter_complex "[0:v]trim=duration=5.375,setpts=PTS-STARTPTS[base];[1:v]scale=720:1280,format=rgba[shade];[base][shade]overlay=0:0:format=auto:shortest=1,tpad=stop_mode=clone:stop_duration=0.625,format=yuv420p[outv];[0:a]atrim=duration=6,asetpts=PTS-STARTPTS,afade=t=out:st=5.65:d=0.35[outa]" \
  -map "[outv]" \
  -map "[outa]" \
  -t 6 \
  -r 24 \
  -c:v libx264 \
  -preset medium \
  -crf 18 \
  -c:a aac \
  -b:a 192k \
  -movflags +faststart \
  "${output_file}"

ffprobe -v error \
  -show_entries format=duration,size:stream=index,codec_type,codec_name,width,height,pix_fmt,r_frame_rate \
  -of default=noprint_wrappers=1 \
  "${output_file}"
