#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
guide_dir="$(cd -- "${script_dir}/.." && pwd)"
drafts_dir="${guide_dir}/renders/flow/drafts"
edits_dir="${guide_dir}/renders/flow/edits"
output_dir="${guide_dir}/renders/edit"
output_file="${1:-${output_dir}/navier-stokes-four-flow-clips-v01-32s.mp4}"

s01="${drafts_dir}/S01_flow_lite_720p_take01.mp4"
s02="${drafts_dir}/S02_flow_lite_720p_take02.mp4"
s05="${drafts_dir}/S05_flow_lite_720p_take01.mp4"
s06="${edits_dir}/S06_flow_lite_take01_8s_label_removed.mp4"

for input_file in "${s01}" "${s02}" "${s05}" "${s06}"; do
  if [[ ! -f "${input_file}" ]]; then
    echo "Missing required input: ${input_file}" >&2
    exit 1
  fi
done

mkdir -p -- "${output_dir}" "$(dirname -- "${output_file}")"

# Four complete Flow generations in story order. No trim, retime, freeze,
# transition, title, placeholder or speed change. Source audio is retained.
# S06 is the full source with only its generated percentage label covered.
ffmpeg -loglevel error -y \
  -i "${s01}" \
  -i "${s02}" \
  -i "${s05}" \
  -i "${s06}" \
  -filter_complex "[0:v]fps=24,scale=720:1280,setsar=1,setpts=PTS-STARTPTS[v0];[0:a]aresample=48000,asetpts=PTS-STARTPTS[a0];[1:v]fps=24,scale=720:1280,setsar=1,setpts=PTS-STARTPTS[v1];[1:a]aresample=48000,asetpts=PTS-STARTPTS[a1];[2:v]fps=24,scale=720:1280,setsar=1,setpts=PTS-STARTPTS[v2];[2:a]aresample=48000,asetpts=PTS-STARTPTS[a2];[3:v]fps=24,scale=720:1280,setsar=1,setpts=PTS-STARTPTS[v3];[3:a]aresample=48000,asetpts=PTS-STARTPTS[a3];[v0][a0][v1][a1][v2][a2][v3][a3]concat=n=4:v=1:a=1[video][audio]" \
  -map "[video]" -map "[audio]" \
  -t 32 \
  -c:v libx264 -preset medium -crf 18 -profile:v high -level 4.0 \
  -c:a aac -b:a 192k -ar 48000 \
  -pix_fmt yuv420p -movflags +faststart \
  "${output_file}"

ffprobe -v error \
  -show_entries format=duration,size:stream=index,codec_type,codec_name,width,height,pix_fmt,r_frame_rate,sample_rate,channels \
  -of default=noprint_wrappers=1 \
  "${output_file}"
