from pathlib import Path
from matplotlib.pylab import f
import numpy as np
from io import StringIO
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec 
from scipy.optimize import curve_fit

plt.rcParams['font.size'] = 24.0
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelsize'] = 'medium'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['lines.linewidth'] = 2.0

def read_block_file(filepath):
    blocks = []

    current_header = []
    current_data = []

    def finalize_block():
        nonlocal current_header, current_data, blocks

        if not current_header and not current_data:
            return

        header_items = []
        for line in current_header:
            line = line[1:].strip()  
            if ":" in line:
                key, value = line.split(":", 1)
                header_items.append((key.strip(), value.strip()))
            else:
                header_items.append((line, None))

        header_dict = dict(header_items)

        if current_data:
            data_str = "\n".join(current_data)
            data = np.loadtxt(StringIO(data_str))
        else:
            data = np.array([])

        blocks.append({
            "header_list": header_items,               "header": header_dict,         
            "data": data
        })

        current_header = []
        current_data = []

    with open(filepath, "r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith("#"):
                if current_data:
                    finalize_block()
                current_header.append(line)
            else:
                current_data.append(line)

    finalize_block()
    return blocks

blocks = read_block_file("C:/VS_Studio_Latex/Uni (lLokale Kopie)/Praktikum/STM_angepasst/Rasterzeit/Au_1s_neueMess.txt")



for i, block in enumerate(blocks):
   for key, value in block["header_list"]:
        print(f"  {key}: {value}")

    #print(block["header"].get("Channel"))

    #print(f"{np.mean(block['data']):.2e} +- {np.std(block['data']):.2e}")
    
N = 187


block_m_1 = np.array(blocks[1]['data'][N]-np.mean(blocks[1]['data'][N]))
block_m_3 = np.array(blocks[3]['data'][N]-np.mean(blocks[3]['data'][N]))
"""
s = np.arange(60, 120)

L1, L1_uns = curve_fit(lambda s, a, b: a * s + b, s, block_m_1[s], p0=[1e-9, 0])
L2, L2_uns = curve_fit(lambda s, a, b: a * s + b, s, block_m_3[s], p0=[1e-9, 0])

f_fit_1 = L1[0] * s + L1[1]
f_fit_2 = L2[0] * s + L2[1]
"""
plt.figure(1, figsize=(20,10), layout='tight')
gs = gridspec.GridSpec(3, 1, height_ratios=[2, 1, 0.1], hspace=0)


plt.subplot(gs[0])
plt.plot(block_m_1, linestyle="--", marker="o", label=blocks[1]['header'].get('Channel'))
plt.plot(block_m_3, linestyle="--", marker="o", label=blocks[3]['header'].get('Channel'))
plt.grid()
#plt.plot(s, f_fit_1, label=f"Flanke vorwärts, a={L1[0]:.2e}±{np.sqrt(np.diag(L1_uns)[0]):.2e}")
#plt.plot(s, f_fit_2, label=f"Flanke rückwärts, a={L2[0]:.2e}±{np.sqrt(np.diag(L2_uns)[0]):.2e}")
plt.ylabel("I/A")
plt.title("Rasterzeit = 1s")
plt.legend()

plt.subplot(gs[1])
plt.grid()
plt.axhline(0, color="black", linestyle="--")
plt.plot((block_m_1-block_m_3), marker="o", linestyle='None')
plt.xlabel("Punkt auf Profillinie")
plt.ylabel(r"$I_{hin}-I_{rück}$")
#plt.show()




#plt.savefig('C:/VS_Studio_Latex/Uni (lLokale Kopie)/Praktikum\STM_angepasst/Rasterzeit/Au_1s_neueMess_I.svg', dpi=300)


M = 1/len(blocks[3]['data']) * np.sum(((blocks[1]['data'][N])-(blocks[3]['data'][N]))**2)

print(f"Mittlere quadratische Abweichung: {M:.2e} m^2")
print(f"Mittelwert I_hin: {np.mean(blocks[0]['data'][N]):.2e} +- {np.std(blocks[0]['data'][N]):.2e} A")
print(f"Mittelwert I_rück: {np.mean(blocks[2]['data'][N]):.2e} +- {np.std(blocks[2]['data'][N]):.2e} A")