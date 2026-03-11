from Gain_auswertung import *
plt.figure(figsize=(10, 8))
farben = ['darkblue', 'darkgreen', 'darkred']
i = 0
for data, U in zip(Daten_namen, [U0, U_80_20, U_60_40]):
    Gain, Gain_uns = plotting(data, U)
    plot_errorfit(data, U, Gain, Gain_uns, c = farben[i])
    plt.yscale('log')
    i+=1
plt.title('Gain bei verschiedenen Gasmischungen')
plt.legend([r"$\pm 1 \sigma$ auf die Anpassung", r"70/30 $CO_2$", r"$\pm 1 \sigma$ auf die Anpassung", r"60/40 $CO_2$", r"$\pm 1 \sigma$ auf die Anpassung", r"80/20 $CO_2$", 'Datenpunkte mit Unsicherheit'], loc = 'lower right')
plt.savefig('Bilder/Verstärkung_alle.svg')
