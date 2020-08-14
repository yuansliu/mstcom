nohup bash -c '/usr/bin/time -v ./minicom -r SRR065389_1.fastq -t 24' </dev/null 1> SRR065389_1.out 2> SRR065389_1.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR445718.fastq -t 24' </dev/null 1> SRR445718.out 2> SRR445718.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR445724.fastq -t 24' </dev/null 1> SRR445724.out 2> SRR445724.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR1265495_1.fastq -t 24' </dev/null 1> SRR1265495_1.out 2> SRR1265495_1.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR1294116.fastq -t 24' </dev/null 1> SRR1294116.out 2> SRR1294116.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR490961.fastq -t 24' </dev/null 1> SRR490961.out 2> SRR490961.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR490976.fastq -t 24' </dev/null 1> SRR490976.out 2> SRR490976.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR554369_1.fastq -t 24 -k 17 -w 3 -m 20' </dev/null 1> SRR554369_1.out 2> SRR554369_1.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR635193_1.fastq -t 24' </dev/null 1> SRR635193_1.out 2> SRR635193_1.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR689233_1.fastq -t 24 -k 17' </dev/null 1> SRR689233_1.out 2> SRR689233_1.err
nohup bash -c '/usr/bin/time -v ./minicom -r ERR532393_1.fastq -t 24' </dev/null 1> ERR532393_1.out 2> ERR532393_1.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR1313062_1.fastq -t 24' </dev/null 1> SRR1313062_1.out 2> SRR1313062_1.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR870667_1.fastq -t 24 -k 25 -m 30 -w 20 -e 18 -S 5' </dev/null 1> SRR870667_1.out 2> SRR870667_1.err
nohup bash -c '/usr/bin/time -v ./minicom -r ERR174310_1.fastq -t 24 -k 25 -m 25 -w 15' </dev/null 1> ERR174310_1.out 2> ERR174310_1.err

nohup bash -c '/usr/bin/time -v ./minicom -1 SRR065389_1.fastq -2 SRR065389_2.fastq -t 24' </dev/null 1> SRR065389_pe.out 2> SRR065389_pe.err
nohup bash -c '/usr/bin/time -v ./minicom -1 ERR532393_1.fastq -2 ERR532393_2.fastq -t 24' </dev/null 1> ERR532393_pe.out 2> ERR532393_pe.err
nohup bash -c '/usr/bin/time -v ./minicom -1 SRR689233_1.fastq -2 SRR689233_2.fastq -t 24 -k 17' </dev/null 1> SRR689233_pe.out 2> SRR689233_pe.err
nohup bash -c '/usr/bin/time -v ./minicom -1 SRR635193_1.fastq -2 SRR635193_2.fastq -t 24' </dev/null 1> SRR635193_pe.out 2> SRR635193_pe.err
nohup bash -c '/usr/bin/time -v ./minicom -1 SRR554369_1.fastq -2 SRR554369_2.fastq -t 24 -k 17 -w 3 -m 20' </dev/null 1> SRR554369_pe.out 2> SRR554369_pe.err
nohup bash -c '/usr/bin/time -v ./minicom -1 SRR1313062_1.fastq -2 SRR1313062_2.fastq -t 24' </dev/null 1> SRR1313062_pe.out 2> SRR1313062_pe.err
nohup bash -c '/usr/bin/time -v ./minicom -1 ERR174310_1.fastq -2 ERR174310_2.fastq -t 24 -t 24 -k 25 -e 4 -w 15 -m 25' </dev/null 1> ERR174310_pe.out 2> ERR174310_pe.err

nohup bash -c '/usr/bin/time -v ./minicom -r SRR065389_1.fastq -p -t 24' </dev/null 1> SRR065389_1.order.out 2> SRR065389_1.order.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR445718.fastq -p -t 24' </dev/null 1> SRR445718.order.out 2> SRR445718.order.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR445724.fastq -p -t 24' </dev/null 1> SRR445724.order.out 2> SRR445724.order.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR1265495_1.fastq -p -t 24' </dev/null 1> SRR1265495_1.order.out 2> SRR1265495_1.order.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR1294116.fastq -p -t 24' </dev/null 1> SRR1294116.order.out 2> SRR1294116.order.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR490961.fastq -p -t 24' </dev/null 1> SRR490961.order.out 2> SRR490961.order.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR490976.fastq -p -t 24' </dev/null 1> SRR490976.order.out 2> SRR490976.order.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR554369_1.fastq -p -t 24 -k 17 -w 3 -m 20' </dev/null 1> SRR554369_1.order.out 2> SRR554369_1.order.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR635193_1.fastq -p -t 24' </dev/null 1> SRR635193_1.order.out 2> SRR635193_1.order.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR689233_1.fastq -p -t 24 -k 17' </dev/null 1> SRR689233_1.order.out 2> SRR689233_1.order.err
nohup bash -c '/usr/bin/time -v ./minicom -r ERR532393_1.fastq -p -t 24' </dev/null 1> ERR532393_1.order.out 2> ERR532393_1.order.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR1313062_1.fastq -p -t 24' </dev/null 1> SRR1313062_1.order.out 2> SRR1313062_1.order.err
nohup bash -c '/usr/bin/time -v ./minicom -r SRR870667_1.fastq -p -t 24 -k 25 -m 30 -w 20 -e 18 -S 5' </dev/null 1> SRR870667_1.order.out 2> SRR870667_1.order.err
nohup bash -c '/usr/bin/time -v ./minicom -r ERR174310_1.fastq -p -t 24 -k 25 -m 25 -w 15' </dev/null 1> ERR174310_1.order.out 2> ERR174310_1.order.err

# decompression
nohup bash -c '/usr/bin/time -v ./minicom -d SRR065389_1_comp.minicom -t 24' </dev/null 1> SRR065389_1_dec.out 2> SRR065389_1_dec.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR445718_comp.minicom -t 24' </dev/null 1> SRR445718_dec.out 2> SRR445718_dec.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR445724_comp.minicom -t 24' </dev/null 1> SRR445724_dec.out 2> SRR445724_dec.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR1265495_1_comp.minicom -t 24' </dev/null 1> SRR1265495_1_dec.out 2> SRR1265495_1_dec.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR1294116_comp.minicom -t 24' </dev/null 1> SRR1294116_dec.out 2> SRR1294116_dec.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR490961_comp.minicom -t 24' </dev/null 1> SRR490961_dec.out 2> SRR490961_dec.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR490976_comp.minicom -t 24' </dev/null 1> SRR490976_dec.out 2> SRR490976_dec.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR554369_1_comp.minicom -t 24' </dev/null 1> SRR554369_1_dec.out 2> SRR554369_1_dec.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR635193_1_comp.minicom -t 24' </dev/null 1> SRR635193_1_dec.out 2> SRR635193_1_dec.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR689233_1_comp.minicom -t 24' </dev/null 1> SRR689233_1_dec.out 2> SRR689233_1_dec.err
nohup bash -c '/usr/bin/time -v ./minicom -d ERR532393_1_comp.minicom -t 24' </dev/null 1> ERR532393_1_dec.out 2> ERR532393_1_dec.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR1313062_1_comp.minicom -t 24' </dev/null 1> SRR1313062_1_dec.out 2> SRR1313062_1_dec.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR870667_1_comp.minicom -t 24' </dev/null 1> SRR870667_1_dec.out 2> SRR870667_1_dec.err
nohup bash -c '/usr/bin/time -v ./minicom -d ERR174310_1_comp.minicom -t 24' </dev/null 1> ERR174310_1_dec.out 2> ERR174310_1_dec.err

nohup bash -c '/usr/bin/time -v ./minicom -d SRR065389_1_comp_order.minicom -t 24' </dev/null 1> SRR065389_1_dec.order.out 2> SRR065389_1_dec.order.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR445718_comp_order.minicom -t 24' </dev/null 1> SRR445718_dec.order.out 2> SRR445718_dec.order.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR445724_comp_order.minicom -t 24' </dev/null 1> SRR445724_dec.order.out 2> SRR445724_dec.order.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR1265495_1_comp_order.minicom -t 24' </dev/null 1> SRR1265495_1_dec.order.out 2> SRR1265495_1_dec.order.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR1294116_comp_order.minicom -t 24' </dev/null 1> SRR1294116_dec.order.out 2> SRR1294116_dec.order.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR490961_comp_order.minicom -t 24' </dev/null 1> SRR490961_dec.order.out 2> SRR490961_dec.order.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR490976_comp_order.minicom -t 24' </dev/null 1> SRR490976_dec.order.out 2> SRR490976_dec.order.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR554369_1_comp_order.minicom -t 24' </dev/null 1> SRR554369_1_dec.order.out 2> SRR554369_1_dec.order.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR635193_1_comp_order.minicom -t 24' </dev/null 1> SRR635193_1_dec.order.out 2> SRR635193_1_dec.order.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR689233_1_comp_order.minicom -t 24' </dev/null 1> SRR689233_1_dec.order.out 2> SRR689233_1_dec.order.err
nohup bash -c '/usr/bin/time -v ./minicom -d ERR532393_1_comp_order.minicom -t 24' </dev/null 1> ERR532393_1_dec.order.out 2> ERR532393_1_dec.order.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR1313062_1_comp_order.minicom -t 24' </dev/null 1> SRR1313062_1_dec.order.out 2> SRR1313062_1_dec.order.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR870667_1_comp_order.minicom -t 24' </dev/null 1> SRR870667_1_dec.order.out 2> SRR870667_1_dec.order.err
nohup bash -c '/usr/bin/time -v ./minicom -d ERR174310_1_comp_order.minicom -t 24' </dev/null 1> ERR174310_1_dec.order.out 2> ERR174310_1_dec.order.err

nohup bash -c '/usr/bin/time -v ./minicom -d SRR065389_comp_pe.minicom -t 24' </dev/null 1> SRR065389_pe_dec.out 2> SRR065389_pe_dec.err
nohup bash -c '/usr/bin/time -v ./minicom -d ERR532393_comp_pe.minicom -t 24' </dev/null 1> ERR532393_pe_dec.out 2> ERR532393_pe_dec.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR689233_comp_pe.minicom -t 24' </dev/null 1> SRR689233_pe_dec.out 2> SRR689233_pe_dec.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR635193_comp_pe.minicom -t 24' </dev/null 1> SRR635193_pe_dec.out 2> SRR635193_pe_dec.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR554369_comp_pe.minicom -t 24' </dev/null 1> SRR554369_pe_dec.out 2> SRR554369_pe_dec.err
nohup bash -c '/usr/bin/time -v ./minicom -d SRR1313062_comp_pe.minicom -t 24' </dev/null 1> SRR1313062_pe_dec.out 2> SRR1313062_pe_dec.err 
nohup bash -c '/usr/bin/time -v ./minicom -d ERR174310_comp_pe.minicom -t 24' </dev/null 1> ERR174310_pe_dec.out 2> ERR174310_pe_dec.err

