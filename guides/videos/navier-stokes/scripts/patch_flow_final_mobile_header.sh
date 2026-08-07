#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
guide_dir="$(cd -- "${script_dir}/.." && pwd)"
input_file="${1:-${guide_dir}/renders/edit/navier-stokes-existence-and-smoothness-final-720p.mp4}"
output_file="${2:-${guide_dir}/renders/edit/navier-stokes-existence-and-smoothness-final-mobile-header.mp4}"
overlay_file="${guide_dir}/04_editorial/overlays/flow-final-header-mobile-patch.ass"

for required_file in "${input_file}" "${overlay_file}"; do
  if [[ ! -f "${required_file}" ]]; then
    echo "Missing required file: ${required_file}" >&2
    exit 1
  fi
done

ffmpeg -hide_banner -y -i "${input_file}" \
  -map 0:v:0 -map 0:a:0 \
  -vf "delogo=x=78:y=105:w=564:h=40:show=0:enable='between(t,7.4,29.7)',delogo=x=58:y=278:w=604:h=38:show=0:enable='between(t,7.4,29.7)',subtitles='${overlay_file}'" \
  -c:v libx264 -preset medium -crf 18 -profile:v high -level 4.0 \
  -c:a copy -movflags +faststart "${output_file}"

ffprobe -v error -count_frames \
  -show_entries format=duration,size:stream=index,codec_name,width,height,r_frame_rate,nb_read_frames,duration,sample_rate,channels \
  -of default=noprint_wrappers=1 "${output_file}"
