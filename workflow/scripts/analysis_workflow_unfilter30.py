
import os
import numpy as np
# import matplotlib.pyplot as plt


def make_demilo_bed(sample,array_file,CENP_A_file,IgG_file,outdir,window_size=5000):
    arrays = {}
    with open(array_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            if sample not in items[0]:
                continue
            print(items[0])
            array = items[0].split(':')
            contig = array[0]
            if contig not in arrays.keys():
                start = int(array[1].split('-')[0])
                end = int(array[1].split('-')[1])
                arrays[contig] = [start,end,items[0]]
                
    CENP_A = {}
    print('reads cenp-a')
    with open(CENP_A_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            start = int(items[1]) 
            if items[3] != 'a':
                continue

            if items[0] not in CENP_A.keys():
                CENP_A[items[0]] = []    
            CENP_A[items[0]].append([start,float(items[10])])
    

    IgG = {}
    print('read IgG')
    with open(IgG_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            start = int(items[1]) 
            if items[3] != 'a':
                continue

            if items[0] not in IgG.keys():
                IgG[items[0]] = []    
            IgG[items[0]].append([start,float(items[10])])

    for i in arrays.keys():
        print(i)
        contig = i
        region = arrays[i]
        size = region[1] - region[0]
        array = region[2]
        print(region)
        numer = int(size / window_size) + 1

        CENP_A_window = {}
        for j in range(numer):
            CENP_A_window[j] = []
        
        contig_CENP_A = CENP_A[i]
        contig_IgG = IgG[i]
        
        for j in contig_CENP_A:
            start = int(j[0]) - region[0]
            if start < 0:
                continue
            index = int(start / window_size)
            value = j[1]
            if start < 0:
                continue
            if index >= len(CENP_A_window.keys()):
                continue
            CENP_A_window[index].append(float(value))
        
        avg_CENPA_window = {}
        for j in CENP_A_window.keys():
            if len(CENP_A_window[j]) == 0:
                avg_CENPA_window[j] = 0
            else:
                avg_CENPA_window[j] = np.mean(CENP_A_window[j])

        
        IgG_window = {}
        for j in range(numer):
            IgG_window[j] = []
        
        for j in contig_IgG:
            start = int(j[0]) - region[0]
            if start < 0:
                continue
            index = int(start / window_size)
            value = j[1]
            if start < 0:
                continue
            if index >= len(IgG_window.keys()):
                continue
            IgG_window[index].append(float(value))
        
        avg_IgG_window = {}
        for j in IgG_window.keys():
            if len(IgG_window[j]) == 0:
                avg_IgG_window[j] = 0
            else:
                avg_IgG_window[j] = np.mean(IgG_window[j])
        
        x = []
        y = []
        index = 0
        max_number = -10000000
        for j in avg_CENPA_window.keys():
            x.append(index)
            y.append(avg_CENPA_window[j] - avg_IgG_window[j])

            index += 1
        
        
        
        outfile = outdir + '/' +  contig + '.xls'
        outfile = open(outfile,'w')
        #outfile.write('index\tcount\n')
        
        for j in range(len(x)):
            start = region[0] + j * window_size
            end = region[0] + j * window_size + window_size
            # outfile.write(array + '\t' + str(start) + '\t' + str(end) + '\t' + str(y[j] / max_number) + '\n')
            outfile.write(array + '\t' + str(start) + '\t' + str(end) + '\t' + str(y[j]) + '\n')
        
        outfile.close()

def read_CENPA(file):
    contig_data = {}
    with open(file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            if items[0] not in contig_data.keys():
                contig_data[items[0]] = []
            contig_data[items[0]].append([int(items[1]),int(items[2]),items[3]])
    return contig_data

def plots(cenplot_dir,DeMiLo_dir,cenpa_file):
    # 获得demi
    cut_and_run_flag = False
    cenpa_data = {}
    if cenpa_file != '':
        cut_and_run_flag = True
        cenpa_data = read_CENPA(cenpa_file)
    
    dirs = os.listdir(cenplot_dir)
    for i in dirs:
        print('plot' + i)
        contig = i.split(':')[0]
        outdir = cenplot_dir + '/' + i
        if cut_and_run_flag == True:
            cenpa = cenpa_data[contig]
            outfile = outdir + '/CutAndRun.bed'
            sorted_cenpa = sorted(cenpa,key=lambda x:x[0])
            outfile = open(outfile,'w')
            array_start = int(i.split(':')[1].split('-')[0])
            array_end = int(i.split(':')[1].split('-')[1])
            max_value = -1
            min_value = 10000000
            for j in sorted_cenpa:
                if j[0] <= array_start:
                    continue
                if j[1] >= array_end:
                    continue
                if float(j[2]) > max_value:
                    max_value = float(j[2])
                if float(j[2]) < min_value:
                    min_value = float(j[2])
            for j in sorted_cenpa:
                if j[0] <= array_start:
                    continue
                if j[1] >= array_end:
                    continue
                # value = float(j[2]) / max_value
                value = float(j[2])
                outfile.write(i + '\t' + str(j[0]) + '\t' + str(j[1]) + '\t' + str(value) + '\n')
            outfile.close()
        
        file = DeMiLo_dir + '/' + contig + '.xls'
        outfile = outdir + '/demilo.bed'
        outfile = open(outfile,'w')
        
        with open(file,'r') as f:
            while True:
                line = f.readline()[:-1]
                if not line:
                    break
                items = line.split('\t')
                if float(items[-1]) < 0:
                    outfile.write(items[0] + '\t' + items[1] + '\t' + items[2] + '\t0\n')
                else:
                    outfile.write(line + '\n')
        
        outfile.close()        
        
        cmd = 'cp ' +  DeMiLo_dir + '/' + contig + '.xls ' + ' ' + outdir + '/demilo.bed'
        os.system(cmd)
        
        # make toml
        tomlfile_name = outdir + '/' + i + '.toml'
        tomlfile = open(tomlfile_name,'w')
        
        # [[tracks]]
        # position = "relative"
        # type = "label"
        # proportion = 0.0025
        # path = "/project/logsdon_shared/projects/HGSVC3/SG_working/demilo/NA20355/NA20355_chr10_haplotype1-0000010:38512490-41842841/cdrs.bed"
        # [tracks.options]
        # color = "black"
        # legend = false
        # hide_x = true
        tomlfile.write('[[tracks]]\n')
        tomlfile.write('position = "relative"\n')
        tomlfile.write('type = "label"\n')
        tomlfile.write('proportion = 0.0025\n')
        tomlfile.write('path = "' + outdir + '/cdrs.bed"\n')
        tomlfile.write('[tracks.options]\n')
        tomlfile.write('color = "black"\n')
        tomlfile.write('legend = false\n')
        tomlfile.write('hide_x = true\n')
        tomlfile.write('\n')
        
        # [[tracks]]
        # position = "relative"
        # proportion = 0.025
        # type = "label"
        # path = "/project/logsdon_shared/projects/HGSVC3/SG_working/demilo/NA20355/NA20355_chr10_haplotype1-0000010:38512490-41842841/sat_annot.bed"
        # [tracks.options]
        # legend = true
        # hide_x = true
        # legend_title = "Structure"
        # legend_title_only = true
        tomlfile.write('[[tracks]]\n')
        tomlfile.write('position = "relative"\n')
        tomlfile.write('proportion = 0.025\n')
        tomlfile.write('type = "label"\n')
        tomlfile.write('path = "' + outdir + '/sat_annot.bed"\n')
        tomlfile.write('[tracks.options]\n')
        tomlfile.write('legend = true\n')
        tomlfile.write('hide_x = true\n')
        tomlfile.write('legend_title = "Structure"\n')
        tomlfile.write('legend_title_only = true\n')
        tomlfile.write('\n')
        
        # [[tracks]]
        # position = "overlap"
        # type = "hor"
        # path = "/project/logsdon_shared/projects/HGSVC3/SG_working/demilo/NA20355/NA20355_chr10_haplotype1-0000010:38512490-41842841/stv.bed"
        # [tracks.options]
        # legend = false
        # sort_order = "descending"
        # hide_x = true
        # bg_border = false
        # mer_filter = 0
        # hor_filter = 0
        # live_only = false
        tomlfile.write('[[tracks]]\n')
        tomlfile.write('position = "overlap"\n')
        tomlfile.write('type = "hor"\n')
        tomlfile.write('path = "' + outdir + '/stv.bed"\n')
        tomlfile.write('[tracks.options]\n')
        tomlfile.write('legend = false\n')
        tomlfile.write('sort_order = "descending"\n')
        tomlfile.write('hide_x = true\n')
        tomlfile.write('bg_border = false\n')
        tomlfile.write('mer_filter = 0\n')
        tomlfile.write('hor_filter = 0\n')
        tomlfile.write('live_only = false\n')
        tomlfile.write('\n')
        
        # [[tracks]]
        # position = "relative"
        # proportion = 0.010
        # type = "strand"
        # path = "/project/logsdon_shared/projects/HGSVC3/SG_working/demilo/NA20355/NA20355_chr10_haplotype1-0000010:38512490-41842841/strand.bed"
        # [tracks.options]
        # legend = false
        # hide_x = true
        # scale = 5
        tomlfile.write('[[tracks]]\n')
        tomlfile.write('position = "relative"\n')
        tomlfile.write('proportion = 0.010\n')
        tomlfile.write('type = "strand"\n')
        tomlfile.write('path = "' + outdir + '/strand.bed"\n')
        tomlfile.write('[tracks.options]\n')
        tomlfile.write('legend = false\n')
        tomlfile.write('hide_x = true\n')
        tomlfile.write('scale = 5\n')
        tomlfile.write('\n')
        
        # [[tracks]]
        # position = "relative"
        # type = "bar"
        # proportion = 0.05
        # path = "/project/logsdon_shared/projects/HGSVC3/SG_working/demilo/NA20355/NA20355_chr10_haplotype1-0000010:38512490-41842841/methyl.bed"
        # [tracks.options]
        # hide_x = true
        # ymax = 1.0
        # ymin = 0
        # legend_title = "Ratio of CpG\nmethylation"
        # legend_title_only = true
        
        tomlfile.write('[[tracks]]\n')
        tomlfile.write('position = "relative"\n')  
        tomlfile.write('type = "bar"\n')
        tomlfile.write('proportion = 0.05\n')
        tomlfile.write('path = "' + outdir + '/methyl.bed"\n')
        tomlfile.write('[tracks.options]\n')
        tomlfile.write('hide_x = true\n')
        tomlfile.write('ymax = 1.0\n')
        tomlfile.write('ymin = 0\n')
        tomlfile.write('legend_title = "Ratio of CpG\\nmethylation"\n')
        tomlfile.write('legend_title_only = true\n')
        tomlfile.write('\n')
        
        # [[tracks]]
        # position = "relative"
        # type = "bar"
        # proportion = 0.05
        # path = "/project/logsdon_shared/projects/HGSVC3/SG_working/demilo/NA20355/NA20355_chr10_haplotype1-0000010:38512490-41842841/demilo.bed"
        # [tracks.options]
        # hide_x = true
        # #ymax = 1.0
        # ymin = 0
        # legend_title = "Normalized CENP-A\nDiMeLo-Seq signals\n(relative to IgG)"
        # legend_title_only = true
        tomlfile.write('[[tracks]]\n')
        tomlfile.write('position = "relative"\n')
        tomlfile.write('type = "bar"\n')
        tomlfile.write('proportion = 0.05\n')
        tomlfile.write('path = "' + outdir + '/demilo.bed"\n')
        tomlfile.write('[tracks.options]\n')
        tomlfile.write('hide_x = true\n')
        tomlfile.write('ymin = 0\n')
        tomlfile.write('legend_title = "Normalized CENP-A\\nDiMeLo-Seq signals\\n(relative to IgG)"\n')
        tomlfile.write('legend_title_only = true\n')
        tomlfile.write('\n')
        
        if cut_and_run_flag == True:
            # [[tracks]]
            # position = "relative"
            # type = "bar"
            # proportion = 0.05
            # path = "/project/logsdon_shared/projects/HGSVC3/SG_working/demilo/NA20355/NA20355_chr10_haplotype1-0000010:38512490-41842841/CutAndRun.bed"
            # [tracks.options]
            # hide_x = true
            # ymin = 0
            # legend_title = "Fold enrichment of\nCENP-A:IgG CUT&RUN data"
            # legend_title_only = true
            tomlfile.write('[[tracks]]\n')
            tomlfile.write('position = "relative"\n')
            tomlfile.write('type = "bar"\n')
            tomlfile.write('proportion = 0.05\n')
            tomlfile.write('path = "' + outdir + '/CutAndRun.bed"\n')
            tomlfile.write('[tracks.options]\n')
            tomlfile.write('hide_x = true\n')
            tomlfile.write('ymin = 0\n')
            tomlfile.write('legend_title = "Fold enrichment of\\nCENP-A:IgG CUT&RUN data"\n')
            tomlfile.write('legend_title_only = true\n')
            tomlfile.write('\n')
        
        # [[tracks]]
        # position = "relative"
        # type = "selfident"
        # proportion = 0.2
        # path = "/project/logsdon_shared/projects/HGSVC3/SG_working/demilo/NA20355/NA20355_chr10_haplotype1-0000010:38512490-41842841/ident.bed"
        # [tracks.options]
        # legend = false
        # hide_x = true
        # invert = true
        tomlfile.write('[[tracks]]\n')
        tomlfile.write('position = "relative"\n')
        tomlfile.write('type = "selfident"\n')
        tomlfile.write('proportion = 0.2\n')
        tomlfile.write('path = "' + outdir + '/ident.bed"\n') 
        tomlfile.write('[tracks.options]\n')
        tomlfile.write('legend = false\n')
        tomlfile.write('hide_x = true\n')
        tomlfile.write('invert = true\n')
        tomlfile.write('\n')
        
        # [[tracks]]
        # position = "relative"
        # type = "position"
        # proportion = 0.005

        # [settings]
        # format = [ "png", "pdf",]
        # transparent = false
        # dim = [ 10, 8,]
        # legend_prop = 0.05
        # axis_h_pad = 0.01
        # dpi = 600
        # legend_pos = "left"
        # layout = "constrained"
        
        tomlfile.write('[[tracks]]\n')
        tomlfile.write('position = "relative"\n')
        tomlfile.write('type = "position"\n')
        tomlfile.write('proportion = 0.005\n')
        tomlfile.write('\n')
        tomlfile.write('[settings]\n')
        tomlfile.write('format = [ "png", "pdf" ]\n')
        tomlfile.write('transparent = false\n')
        tomlfile.write('dim = [ 10, 8 ]\n')
        tomlfile.write('legend_prop = 0.05\n')
        tomlfile.write('axis_h_pad = 0.01\n')
        tomlfile.write('dpi = 600\n')
        tomlfile.write('legend_pos = "left"\n')
        tomlfile.write('layout = "constrained"\n')
        tomlfile.close()
    
    # run
    for i in dirs:
        # /project/logsdon_shared/projects/HGSVC3/SG_working/demilo/other/HG02769_unfilter30/
        # cenplot/HG02769_chr20_haplotype1-0000051:25820384-30734683/HG02769_chr20_haplotype1-0000051:25820384-30734683.png
        outdir = cenplot_dir + '/' + i
        if os.path.exists(outdir + '/' + i + '.png'):
            print(i)
            continue
        cmd = 'cd ' + '/project/logsdon_shared/projects/Keith/CenPlot'
        os.system(cmd)
        print('run plot')
        print(i)
        tomlfile_name = outdir + '/' + i + '.toml'
        # /project/logsdon_shared/projects/Keith/CenPlot/venv/bin/python 
        # -m cenplot.main draw 
        # -t /project/logsdon_shared/projects/Keith/CenPlot/HG00096_chr1_haplotype1-0000018:120846650-126839919.toml 
        # -d . -c HG00096_chr1_haplotype1-0000018:120846650-126839919
        cmd = '/project/logsdon_shared/projects/Keith/CenPlot/venv/bin/python ' + \
            ' -m ' + ' cenplot.main draw ' + ' -t ' + tomlfile_name + ' -d ' + outdir + ' -c ' + i
        print(cmd)
        os.system(cmd)


def main():
    sample = 'HG03065'
    
    CENP_A_file = '/project/logsdon_shared/projects/HGSVC3/SG_working/demilo/other/'+sample+'_unfilter30/'+sample+'_treatment_2.0.bed'
    IgG_file = '/project/logsdon_shared/projects/HGSVC3/SG_working/demilo/other/'+sample+'_unfilter30/'+sample+'_control_2.0.bed'
    cenplot_dir = '/project/logsdon_shared/projects/HGSVC3/SG_working/demilo/other/'+sample+'_unfilter30/cenplot'
    outDeMiLodir = '/project/logsdon_shared/projects/HGSVC3/SG_working/demilo/other/' + sample + '_unfilter30'
    
    cenpa_file = ''
    
    
    if not os.path.exists(outDeMiLodir):
        os.makedirs(outDeMiLodir)
    
    array_file = '/project/logsdon_shared/projects/HGSVC3/SG_working/demilo/all_array_length.bed'
    print('make demilo data')
    make_demilo_bed(sample,array_file,CENP_A_file,IgG_file,outDeMiLodir)
    print('plot')
    plots(cenplot_dir,outDeMiLodir,cenpa_file)

    
        
if __name__ == '__main__':
    main()
