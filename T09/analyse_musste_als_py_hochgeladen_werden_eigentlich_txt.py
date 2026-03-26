//
//  define variables
//

// dummy variables
int n;
int k;

int result;
int ktot; 
float px;
float py;
float pz;
float px2;
float py2;
float pz2;
float masse;

int irecl, istat;
int icycle;

int nevmax = 99999; // max number events
int nevent = 0; // event counter
float s_H = 91.2;
float s_89 = 89.48;
float s_91 = 91.33;
float s_93 = 93.02;
float s = s_93;
// input data file
char datfile[] = "../datasets/profiles/zblc/93gev.dat";
//char datfile[] = "../datasets/profiles/zblc/muons.dat";
//char datfile[] = "../datasets/profiles/zblc/hadrons.dat";


// histogram output file name
char hisfile[] = "Output/93gev.root";
//char hisfile[] = "Output/muons.root";
//char hisfile[] = "Output/hadrons.root";



// open root file in which histograms can be stored
TFile *histofile = new TFile(hisfile,"RECREATE");

// create new histogram
TH1F *histo1 = new TH1F("histo1",
"Histogramm der Winkelverteilungen der Events (93 GeV); cos(#theta); Wahrscheinlichkeit", 50,0.,1.005);

//Histogramm der Einteilchenenergien der Teilchen mit Myonenmasse ; Teilchenenergie normiert auf die halbe Schwerpunktsenergie; 
// Histogramm der Winkelverteilungen der Events (91 GeV); cos(#theta) ;
// Histogramm der Winkelverteilungen der Events (91 GeV); cos(#theta) ;
// normierte Gesamtenergie 
// Histogramm der Gesamtenergie 
// Histogramm der Gesamtenergie mit cut 
// Histogramm der Gesamtenergie mit beiden cuts 
// HELP: 
//   Here the histogram histo1 has been reserved,
//      with 25 bins between -50. and 50.:
//        ->  histo1 is the name of the histogram 
//        -> `x Kompo ....' is the histogram title
//        ->  25 is the number of bins
//        -> -50. is the lower limit (here momentum in GeV) 
//                   of the horizontal axis
//        -> +50. ... ... upper ....................................
//
//    More histograms can be defined as desired.


//
//  loop over all events and do analysis
//
for(n=1; n<=nevmax; n++) {
    
    //  create an event of type cevent 
    cevent event;
    float Ges_Energie = 0;
    float Teilchen_Energie = 0;
    int N_cut = 23;
    float cos_winkel;
    float E_cut_lower = 0.85;
    float E_cut_upper = 1.05;
    int n_count = 0;
    int iz1 = 0;
    int iz2 = 0;
    
    //  read in one event from data file
    result = read_event(datfile, event);

    //  note: read_event returns
    //      0  if ....  (see function definition)
    //     -2 if end of file reached

    // analyse only events which are not "empty"
    if(result==0) {
        //  now loop over all PARTICLES within the current event
        
        // get the number of particles in the event
        ktot = event.number_particles();
        if (n==1) {
            cout << endl << "Event " << n << ":" << endl;
            event.print();
        }
            for(k=1; k<=ktot; k++) {
                px = event.momentum(k,1);
                py = event.momentum(k,2);
                pz = event.momentum(k,3);
                masse = event.mass(k); 
                Teilchen_Energie =sqrt(px*px+py*py+pz*pz+masse*masse);
                Ges_Energie = Teilchen_Energie + Ges_Energie;
                //  do some analysis here
                //  for example select muons like so:
                //
                if(fabs(event.mass(k) - 0.106) < 0.001 && E_cut_upper>Teilchen_Energie/(0.5*s) && Teilchen_Energie/(0.5*s)>E_cut_lower) { 
                      if(iz1 == 0){
                        iz1 = k;
                        n_count++;
                    }
                    else if(iz1 != 0){
                        iz2 = k;
                        n_count++;
                    }
                 }
            //if(ktot>=N_cut){
              //  Ges_Energie = Ges_Energie/s;//Hier ist die Energie!
                //if (Ges_Energie >=E_cut_lower){
                  //  nevent++;
               // }
            }
            if(n_count == 2){
                px = event.momentum(iz1,1);
                py = event.momentum(iz1,2);
                pz = event.momentum(iz1,3);
                masse = event.mass(iz1); 
                Teilchen_Energie =sqrt(px*px+py*py+pz*pz+masse*masse);
                //histo1->Fill(Teilchen_Energie/(0.5*s));
                
                px2 = event.momentum(iz2,1);
                py2 = event.momentum(iz2,2);
                pz2 = event.momentum(iz2,3);
                masse = event.mass(iz2); 
                Teilchen_Energie =sqrt(px2*px2+py2*py2+pz2*pz2+masse*masse);
                //histo1->Fill(Teilchen_Energie/(0.5*s));
                
                cos_winkel = fabs(px*px2+py*py2+pz*pz2)/sqrt((px*px+py*py+pz*pz)*(px2*px2+py2*py2+pz2*pz2));
                nevent++;
                histo1->Fill(cos_winkel);
                
            }   
            
        }
    else if(result==-2) {
        break;
    }
}

//  print out some statistics:
cout << "Number of analysed events: " << nevent << endl;
// save histograms to the file defined above
histofile->Write();