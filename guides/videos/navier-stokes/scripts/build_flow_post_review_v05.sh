#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
guide_dir="$(cd -- "${script_dir}/.." && pwd)"

default_input="/home/raxim/Downloads/Initial_Scene_-_2026-08-05_202608051411.mp4"
input_file="${1:-${default_input}}"
output_file="${2:-${guide_dir}/renders/edit/navier-stokes-flow-post-v05-review-720p.mp4}"
overlay_file="${guide_dir}/04_editorial/overlays/flow-final-overlays-v01.ass"

for required_file in "${input_file}" "${overlay_file}"; do
  if [[ ! -f "${required_file}" ]]; then
    echo "Missing required file: ${required_file}" >&2
    exit 1
  fi
done

mkdir -p -- "$(dirname -- "${output_file}")"

# Keep the V04 visual overlap repair. The left-hand perimeter glint at
# 7.85-8.03 s has been removed from the subtitle overlay.
#
# Audio changes relative to V04:
# - preserve the original audio continuously through 7 s (the old 12 ms
#   fade-out/fade-in created an unnecessary small hesitation);
# - repair only the real discontinuity at exactly 8.000 s with a 40 ms
#   equal-power crossfade, removing that tiny restart fragment;
# - keep the already-clean 15 s and 23 s continuation joins unchanged;
# - compensate the 40 ms splice by an inaudible 0.135% tempo correction so
#   the soundtrack remains exactly aligned to the 29.7 s video.
ffmpeg -hide_banner -y \
  -i "${input_file}" \
  -filter_complex "\
[0:v]split=3[v0][v1][v2];\
[v0]trim=start=0:end=15.000,setpts=PTS-STARTPTS[va];\
[v1]trim=start=15.958333:end=23.000,setpts=(PTS-STARTPTS)*1.136095[vb];\
[v2]trim=start=23.958333:end=31.000,setpts=(PTS-STARTPTS)*1.136095[vc];\
[va][vb][vc]concat=n=3:v=1:a=0,\
fps=30,scale=720:1280:flags=lanczos,setsar=1,\
tpad=stop_mode=clone:stop_duration=0.034,trim=duration=29.700,\
subtitles='${overlay_file}',format=yuv420p[vout];\
[0:a]aresample=48000,asplit=3[a0][a1][a2];\
[a0]atrim=start=0:end=15.000,asetpts=PTS-STARTPTS[pre];\
[a1]atrim=start=15.083333:end=23.000,asetpts=PTS-STARTPTS[p1];\
[pre][p1]acrossfade=d=0.200000:c1=qsin:c2=qsin[mid];\
[a2]atrim=start=23.083333:end=30.150000,asetpts=PTS-STARTPTS[p2];\
[mid][p2]acrossfade=d=0.083333:c1=qsin:c2=qsin,\
apad=pad_dur=0.050,atrim=duration=29.700[abase];\
[abase]asplit=2[aleft][aright];\
[aleft]atrim=start=0:end=8.000,asetpts=PTS-STARTPTS[al];\
[aright]atrim=start=8.000,asetpts=PTS-STARTPTS[ar];\
[al][ar]acrossfade=d=0.040000:c1=qsin:c2=qsin,\
atempo=0.998653,apad=pad_dur=0.050,atrim=duration=29.700,\
volume=-1.2dB[aout]" \
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
