#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
guide_dir="$(cd -- "${script_dir}/.." && pwd)"

default_input="/home/raxim/Downloads/Initial_Scene_-_2026-08-05_202608051411.mp4"
input_file="${1:-${default_input}}"
output_file="${2:-${guide_dir}/renders/edit/navier-stokes-flow-post-v02-review-720p.mp4}"
overlay_file="${guide_dir}/04_editorial/overlays/flow-final-overlays-v01.ass"

for required_file in "${input_file}" "${overlay_file}"; do
  if [[ ! -f "${required_file}" ]]; then
    echo "Missing required file: ${required_file}" >&2
    exit 1
  fi
done

mkdir -p -- "$(dirname -- "${output_file}")"

# Review at the source's native 720 x 1280 resolution. The two persistent
# geometry jumps begin at source frames 360 (15.000 s) and 552 (23.000 s).
# At each boundary drop only the first two frames of the new state and use a
# two-frame dissolve. Do not run a broadband de-click filter: preserve the
# generated music and repair joins only with short fades/crossfades.
ffmpeg -hide_banner -y \
  -i "${input_file}" \
  -filter_complex "\
[0:v]split=3[v0][v1][v2];\
[v0]trim=start=0:end=15.000,setpts=PTS-STARTPTS[va];\
[v1]trim=start=15.083333:end=23.000,setpts=PTS-STARTPTS[vb];\
[va][vb]xfade=transition=fade:duration=0.083333:offset=14.916667[vab];\
[v2]trim=start=23.083333:end=30.033333,setpts=PTS-STARTPTS[vc];\
[vab][vc]xfade=transition=fade:duration=0.083333:offset=22.750000,\
fps=30,scale=720:1280:flags=lanczos,setsar=1,\
tpad=stop_mode=clone:stop_duration=0.034,trim=duration=29.700,\
subtitles='${overlay_file}',format=yuv420p[vout];\
[0:a]aresample=48000,asplit=4[a0][a1][a2][a3];\
[a0]atrim=start=0:end=7.000,asetpts=PTS-STARTPTS,\
afade=t=out:st=6.988:d=0.012[p0];\
[a1]atrim=start=7.000:end=15.000,asetpts=PTS-STARTPTS,\
afade=t=in:st=0:d=0.012[p1];\
[p0][p1]concat=n=2:v=0:a=1[pre];\
[a2]atrim=start=15.083333:end=23.000,asetpts=PTS-STARTPTS[p2];\
[pre][p2]acrossfade=d=0.083333:c1=qsin:c2=qsin[mid];\
[a3]atrim=start=23.083333:end=30.033333,asetpts=PTS-STARTPTS[p3];\
[mid][p3]acrossfade=d=0.083333:c1=qsin:c2=qsin,\
apad=pad_dur=0.050,atrim=duration=29.700,volume=-1.2dB[aout]" \
  -map "[vout]" -map "[aout]" \
  -frames:v 891 -t 29.700 \
  -c:v libx264 -preset medium -crf 18 -profile:v high -level 4.0 \
  -c:a aac -b:a 256k -ar 48000 \
  -movflags +faststart \
  "${output_file}"

ffprobe -v error -count_frames \
  -show_entries format=duration,size:stream=index,codec_name,width,height,r_frame_rate,nb_read_frames,sample_rate,channels \
  -of default=noprint_wrappers=1 \
  "${output_file}"
