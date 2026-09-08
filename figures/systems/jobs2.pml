@render.pml
set_color cpep, [74, 58, 167]
delete all
load /Users/bright/Documents/MOSAICS/work/preprint-I/systems/structures/peptide.pdb, m
hide everything; show sticks, not hydro; color cpep, elem C; set stick_radius, 0.25; show spheres, elem N+O; set sphere_scale, 0.22; color tv_blue, elem N; color tv_red, elem O; hide spheres, hydro
orient m
zoom m, 2
png /private/tmp/claude-501/-Users-bright-Library-Application-Support-Claude-scratch-workspaces-2f1eee86-9420-4fdb-ad9a-58b855c266ac-bee9ec2c-588a-47bb-9545-cfef81db5b89-scratch-2026-09-08-a08093/a173a970-b61e-4d20-bf9e-40c516c0299d/scratchpad/paper/figures/systems/peptide.png, width=2000, height=1500, dpi=300, ray=1
delete all
load /Users/bright/Documents/MOSAICS/work/preprint-I/systems/structures/tax9.pdb, m
hide everything; show sticks, not hydro; color cpep, elem C; set stick_radius, 0.25; show spheres, elem N+O; set sphere_scale, 0.22; color tv_blue, elem N; color tv_red, elem O; hide spheres, hydro
orient m
zoom m, 2
png /private/tmp/claude-501/-Users-bright-Library-Application-Support-Claude-scratch-workspaces-2f1eee86-9420-4fdb-ad9a-58b855c266ac-bee9ec2c-588a-47bb-9545-cfef81db5b89-scratch-2026-09-08-a08093/a173a970-b61e-4d20-bf9e-40c516c0299d/scratchpad/paper/figures/systems/tax9.png, width=2000, height=1500, dpi=300, ray=1
delete all
load /Users/bright/Documents/MOSAICS/work/preprint-I/systems/structures/1HHK.pdb, m
hide everything; show cartoon; color cprot, chain A; color cprot2, chain B; show sticks, chain C and not hydro; color cpep, chain C and elem C; set stick_radius, 0.3
orient m
zoom m, 2
png /private/tmp/claude-501/-Users-bright-Library-Application-Support-Claude-scratch-workspaces-2f1eee86-9420-4fdb-ad9a-58b855c266ac-bee9ec2c-588a-47bb-9545-cfef81db5b89-scratch-2026-09-08-a08093/a173a970-b61e-4d20-bf9e-40c516c0299d/scratchpad/paper/figures/systems/1HHK.png, width=2000, height=1500, dpi=300, ray=1
delete all
load /Users/bright/Documents/MOSAICS/work/preprint-I/systems/structures/5PTI.pdb, m
hide everything; show cartoon; color cprot; show sticks, resn CYX and not hydro; color cprot, resn CYX and elem C; color yellow, resn CYX and elem S; set stick_radius, 0.3
orient m
zoom m, 2
png /private/tmp/claude-501/-Users-bright-Library-Application-Support-Claude-scratch-workspaces-2f1eee86-9420-4fdb-ad9a-58b855c266ac-bee9ec2c-588a-47bb-9545-cfef81db5b89-scratch-2026-09-08-a08093/a173a970-b61e-4d20-bf9e-40c516c0299d/scratchpad/paper/figures/systems/5PTI.png, width=2000, height=1500, dpi=300, ray=1
quit
