#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
video_dir="$(cd -- "$script_dir/.." && pwd)"
frame_dir="$video_dir/assets/previs-10"
output="$video_dir/renders/edit/navier-stokes-10-frame-morph-previs-v01-30s.mp4"

frames=(
  "$frame_dir/K01-00s-opening-world.png"
  "$frame_dir/K02-03s-opening-fast.png"
  "$frame_dir/K03-06s-opening-decelerates.png"
  "$frame_dir/K04-09s-close-cube.png"
  "$frame_dir/K05-12s-glass-boundary.png"
  "$frame_dir/K06-15s-inside-macro.png"
  "$frame_dir/K07-18s-slow-macro-track.png"
  "$frame_dir/K08-21s-deep-test-macro.png"
  "$frame_dir/K09-24s-frontier-heightfield.png"
  "$frame_dir/K10-27s-impact-applications.png"
)

for frame in "${frames[@]}"; do
  if [[ ! -f "$frame" ]]; then
    echo "Missing required previs frame: $frame" >&2
    exit 1
  fi
done

mkdir -p "$(dirname -- "$output")"

ffmpeg -hide_banner -y \
  -loop 1 -t 3 -i "${frames[0]}" \
  -loop 1 -t 3 -i "${frames[1]}" \
  -loop 1 -t 3 -i "${frames[2]}" \
  -loop 1 -t 3 -i "${frames[3]}" \
  -loop 1 -t 3 -i "${frames[4]}" \
  -loop 1 -t 3 -i "${frames[5]}" \
  -loop 1 -t 3 -i "${frames[6]}" \
  -loop 1 -t 3 -i "${frames[7]}" \
  -loop 1 -t 3 -i "${frames[8]}" \
  -loop 1 -t 3 -i "${frames[9]}" \
  -filter_complex \
  "[0:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS[k0];
   [1:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS,split=2[k1a][k1b];
   [2:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS,split=2[k2a][k2b];
   [3:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS,split=2[k3a][k3b];
   [4:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS,split=2[k4a][k4b];
   [5:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS,split=2[k5a][k5b];
   [6:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS,split=2[k6a][k6b];
   [7:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS,split=2[k7a][k7b];
   [8:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS,split=2[k8a][k8b];
   [9:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS,split=2[k9a][k9hold];
   [k0][k1a]xfade=transition=fade:duration=3:offset=0[t01];
   [k1b][k2a]xfade=transition=fade:duration=3:offset=0[t02];
   [k2b][k3a]xfade=transition=fade:duration=3:offset=0[t03];
   [k3b][k4a]xfade=transition=fade:duration=3:offset=0[t04];
   [k4b][k5a]xfade=transition=fade:duration=3:offset=0[t05];
   [k5b][k6a]xfade=transition=fade:duration=3:offset=0[t06];
   [k6b][k7a]xfade=transition=fade:duration=3:offset=0[t07];
   [k7b][k8a]xfade=transition=fade:duration=3:offset=0[t08];
   [k8b][k9a]xfade=transition=fade:duration=3:offset=0[t09];
   [k9hold]trim=duration=3,setpts=PTS-STARTPTS[t10];
   [t01][t02][t03][t04][t05][t06][t07][t08][t09][t10]concat=n=10:v=1:a=0[v]" \
  -map "[v]" -t 30 \
  -c:v libx264 -preset medium -crf 18 -pix_fmt yuv420p \
  -movflags +faststart \
  "$output"

echo "Created: $output"
