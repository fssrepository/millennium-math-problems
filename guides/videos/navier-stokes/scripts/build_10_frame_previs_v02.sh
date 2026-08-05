#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
video_dir="$(cd -- "$script_dir/.." && pwd)"
frame_dir="$video_dir/assets/previs-10"
v02_dir="$video_dir/assets/previs-10-v02"
output="$video_dir/renders/edit/navier-stokes-10-frame-story-previs-v02-30s.mp4"
contact_sheet="$v02_dir/contact-sheet.png"
render_dir="$(mktemp -d)"
trap 'rm -rf -- "$render_dir"' EXIT

# K04–K09 share one source and progressively tighter crops so the planning
# preview reads as one forward camera move. The overlays change; the underlying
# fluid world does not jump to a different still at the far-shell story beat.
k04_raw="$render_dir/K04-09s-continuous-fluid-raw.png"
k05_raw="$render_dir/K05-12s-continuous-fluid-raw.png"
k04_base="$render_dir/K04-09s-continuous-fluid.png"
k05_base="$render_dir/K05-12s-continuous-fluid.png"
pre_application_mask="$render_dir/pre-application-focus-mask.png"
pre_application_veil="$render_dir/pre-application-focus-veil.png"
convert "$frame_dir/K10-27s-impact-applications.png" \
  -crop '400x711+300+512' +repage -resize '720x1280!' "$k04_raw"
convert "$frame_dir/K10-27s-impact-applications.png" \
  -crop '365x649+300+591' +repage -resize '720x1280!' "$k05_raw"
convert -size 720x1280 xc:'gray97' -fill black \
  -draw 'roundrectangle -40,390 630,1100 60,60' -blur 0x45 \
  "$pre_application_mask"
convert -size 720x1280 xc:'#02070d' "$pre_application_mask" -alpha off \
  -compose copy_opacity -composite "$pre_application_veil"
convert "$k04_raw" "$pre_application_veil" -compose over -composite "$k04_base"
convert "$k05_raw" "$pre_application_veil" -compose over -composite "$k05_base"

for source in "$v02_dir"/K{06,07,08,09,10}-*.svg; do
  overlay="$render_dir/$(basename -- "${source%.svg}")-overlay.png"
  prepared_base="$render_dir/$(basename -- "${source%.svg}")-base.png"
  rendered="$render_dir/$(basename -- "${source%.svg}").png"
  (
    cd -- "$v02_dir"
    convert -background none "$(basename -- "$source")" "$overlay"
  )
  case "$(basename -- "$source")" in
    K06-*) crop='336x597+302+650' ;;
    K07-*) crop='277x492+332+710' ;;
    K08-*) crop='220x391+360+760' ;;
    K09-*) crop='170x302+385+805' ;;
    K10-*) crop='' ;;
  esac
  if [[ -n "$crop" ]]; then
    convert "$frame_dir/K10-27s-impact-applications.png" -crop "$crop" +repage \
      -resize '720x1280!' "$prepared_base"
  else
    convert "$frame_dir/K10-27s-impact-applications.png" \
      -resize '720x1280^' -gravity center -extent 720x1280 "$prepared_base"
  fi
  convert "$prepared_base" "$overlay" \
    -compose over -composite "$rendered"
done

frames=(
  "$frame_dir/K01-00s-opening-world.png"
  "$frame_dir/K02-03s-opening-fast.png"
  "$frame_dir/K03-06s-opening-decelerates.png"
  "$k04_base"
  "$k05_base"
  "$render_dir/K06-15s-expected-decay.png"
  "$render_dir/K07-18s-persistent-response.png"
  "$render_dir/K08-21s-obstruction.png"
  "$render_dir/K09-24s-save-refine.png"
  "$render_dir/K10-27s-continuous-impact.png"
)

for frame in "${frames[@]}"; do
  if [[ ! -f "$frame" ]]; then
    echo "Missing required v02 previs frame: $frame" >&2
    exit 1
  fi
done

mkdir -p "$(dirname -- "$output")"

montage \
  "${frames[0]}" "${frames[1]}" "${frames[2]}" "${frames[3]}" "${frames[4]}" \
  "${frames[5]}" "${frames[6]}" "${frames[7]}" "${frames[8]}" "${frames[9]}" \
  -thumbnail '216x384^' -gravity center -extent 216x384 \
  -background '#05090d' -gravity center -tile 5x2 -geometry +12+12 "$contact_sheet"

ffmpeg -hide_banner -y \
  -loop 1 -t 3.5 -i "${frames[0]}" \
  -loop 1 -t 3.5 -i "${frames[1]}" \
  -loop 1 -t 3.5 -i "${frames[2]}" \
  -loop 1 -t 3.5 -i "${frames[3]}" \
  -loop 1 -t 3.5 -i "${frames[4]}" \
  -loop 1 -t 3.5 -i "${frames[5]}" \
  -loop 1 -t 3.5 -i "${frames[6]}" \
  -loop 1 -t 3.5 -i "${frames[7]}" \
  -loop 1 -t 3.5 -i "${frames[8]}" \
  -loop 1 -t 3.5 -i "${frames[9]}" \
  -filter_complex \
  "[0:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS[k0];
   [1:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS[k1];
   [2:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS[k2];
   [3:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS[k3];
   [4:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS[k4];
   [5:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS[k5];
   [6:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS[k6];
   [7:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS[k7];
   [8:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS[k8];
   [9:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS,
        zoompan=z='max(1,5.535-on*0.06387)':x='max(0,min(iw-iw/zoom,360-iw/(2*zoom)))':y='max(0,min(ih-ih/zoom,732-ih/(2*zoom)))':d=1:s=720x1280:fps=24,
        trim=duration=3,setpts=PTS-STARTPTS[k9zoom];
   [k0][k1]xfade=transition=fade:duration=0.5:offset=3[t01];
   [t01][k2]xfade=transition=fade:duration=0.5:offset=6[t02];
   [t02][k3]xfade=transition=fade:duration=0.5:offset=9[t03];
   [t03][k4]xfade=transition=fade:duration=0.5:offset=12[t04];
   [t04][k5]xfade=transition=fade:duration=0.5:offset=15[t05];
   [t05][k6]xfade=transition=fade:duration=0.5:offset=18[t06];
   [t06][k7]xfade=transition=fade:duration=0.5:offset=21[t07];
   [t07][k8]xfade=transition=fade:duration=0.5:offset=24[t08];
   [t08]trim=duration=27,setpts=PTS-STARTPTS[main];
   [main][k9zoom]concat=n=2:v=1:a=0[v]" \
  -map "[v]" -t 30 \
  -c:v libx264 -preset medium -crf 18 -pix_fmt yuv420p \
  -movflags +faststart \
  "$output"

echo "Created: $output"
echo "Created: $contact_sheet"
