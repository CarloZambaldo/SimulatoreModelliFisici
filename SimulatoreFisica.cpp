 *                                              *
 *      PROGRAMMA PER SIMULAZIONI FISICHE       *
 *                                              *
 *                written by:                   *
 *              CARLO ZAMBALDO                  *
 *         (carlo.zambaldo@gmail.com)           *
 *                                              *
 *      LICEO GIROLAMO FRACASTORO, VERONA       *
 *                                              *
 *       Last Update: NOVEMBER 06, 2018         *
 *                                              *
 * -------------------------------------------- */

#include <math.h>
#include <fstream>
#include <iostream>
using namespace std;

#define G 6.67408e-11		// costante gravitazionale
#define Mt 5.972e24			// massa della Terra [kg]
#define Rt 6.371e6 			// raggio della Terra [m]
#define epsi 8.854188e-12	// costante dielettrica nel vuoto
#define pi 3.1415926		// valore di pi greco
double fq=1.00;				// frequenza di registrazione dei dati [s]

#define versione "3.0.7"	// Versione del programma

// -------------------------------------------------------------

void prosDat(int scelta); // vedi giu'

char NomeFile(){ // questa funzione dovrebbe modificare il nome del file per evitare omonimie (NON FUNZIONA)
	char Num;
	cout<<"Inserire il numero della simulazione: ";
	cin>>Num;
	return Num;
}
double Freq(double fre, double v){	// Algoritmo per aumentare la frequenza dei risultati quando la velocita' dell'oggetto aumenta (NON FUNZIONA)
	fre=(fre/v)*100;
	return fre;
}

// -------------------------------------------------------------


// MOTO PARABOLICO
void moti(double x, double y, double v, double teta){
	double vx = v*cos(teta), vy = v*sin(teta), h=y, yp=y, j, time; // teta e' in rad, quando viene stampato e' in gradi: (teta*180.0/pi)
	double vyp=vy, g=-(G*Mt)/pow((Rt+y),2.0);
	ofstream Ris1;
	Ris1.open("MotoParabolico.txt", ios::out);
	Ris1<<"Dati inseriti: x = "<<x<<"; y = "<<y<<"; v = "<<v<<"; teta = "<<teta<<endl;
	
	Ris1<<"T [s]\tx [m]\ty [m]\tv [m/s]\tangle\tvx [m/s]\tvy [m/s]\tg [m/s^2]\t"<<endl;
	Ris1<<"0\t"<<x<<"\t"<<y<<"\t"<<v<<"\t"<<(teta*180.0/pi)<<"\t"<<vx<<"\t"<<vy<<"\t"<<g<<"\t"<<endl;
	
	for(double t=1; y>=0; t+=fq){
		// calcolo la velocità di fuga dell'oggetto (se velocità iniziale è maggiore allora break;
		if(vy<=sqrt((2*G*Mt)/(Rt+y))){
			// calcolo a t+=fq ogni istante di tempo varie cose.
			g=-(G*Mt)/pow((Rt+y),2.0);
			x=vx*t;
			y=yp+vyp*t+0.5*g*pow(t,2.0);
			vy=vyp+g*t;
			v=sqrt(pow(vx,2.0)+pow(vy,2.0));
			teta=tan(vy/vx);
			Ris1<<t<<"\t"<<x<<"\t"<<y<<"\t"<<v<<"\t"<<(teta*180.0/pi)<<"\t"<<vx<<"\t"<<vy<<"\t"<<g<<"\t"<<endl; // Stampa su file
		}else{
			cout<<"POSIZIONE FINALE DELL'OGGETTO"<<endl;
			cout<<"L'oggetto, trascurando gli attriti, non arrivera' mai a terra."<<endl;
			Ris1<<"L'oggetto, trascurando gli attriti, non arrivera' mai a terra."<<endl;
			Ris1.close();
			system("PAUSE");
			break;
		}
	}
	
			cout<<"File salvato con successo."<<endl;
			Ris1<<"------------------ END OF FILE ------------------"<<endl;
			Ris1.close();
		
			cout<<"POSIZIONE FINALE DELL'OGGETTO"<<endl;
			// Calcolo la posizione finale dell'oggetto
			
			// bug fixed: j scambiata con tempo, nuova variabile time
			y=0;
			time=sqrt(-2*h/g);
			g=-(G*Mt)/pow((Rt+y),2.0);
			j=(-vyp+sqrt(pow(vyp,2.0)-2*g*h))/(-g);
			x=vx*time;
			vy=vyp+g*time;
			v=sqrt(pow(vx,2.0)+pow(vy,2.0));
			teta=tan(vy/vx);
			
			// Stampo a video la posizione finale dell'oggetto
			cout<<"T [s]\tx [m]\ty [m]\tv [m/s]\tangle"<<endl;
			cout<<time<<"\t"<<x<<"\t"<<y<<"\t"<<v<<"\t"<<(teta*180.0/pi)<<endl;
			system("PAUSE");
}

// FORZA GRAVITAZIONALE - con moto associato
void grav(double A, double B, double d){
	double vA=0, vB=0, pA=A*vA, pB=B*vB, F=(G*A*B)/pow(d,2.0), h=0;
	double fq=1.00;
		ofstream Ris2;
		Ris2.open("ForzaGravitazionale.txt", ios::out);

		Ris2<<"T [s]\td [m]\tp(A) [kg*m/s]\tp(B) [kg*m/s]\tF [N]"<<endl;

		for(double t=0; d>=0; t+=fq){
			Ris2<<t<<"\t"<<d<<"\t"<<pA<<"\t"<<pB<<"\t"<<F<<"\tva:"<<vA<<"\tvb:"<<vB<<endl; // Stampa su file
			F=((G*A*B)/pow(d,2.0));
			vA=(F/A)*t;   pA=A*vA;
			vB=(F/B)*t;   pB=B*vB;

			h=(G*A*B)/((0.5*A*pow(vA,2.0))+(0.5*B*pow(vB,2.0))-((G*A*B)/d))+d;
			d+=h;
		}

		// Da inserire in futuro il richiamo alla procedura degli urti: quando i due corpi sono a distanza d=0 (oppure d=r1+r2)
		// Da richiedere prima se l'urto sara' elastico o meno, poi inviare alla procedura corrispondente
		{
			char ooo;
			cout<<"Vuoi proseguire con l'urto delle due particelle? [S]i - [N]o.";
			if((ooo == 'S')||(ooo == 's')) prosDat(2);
		}
		
		cout<<"File salvato con successo."<<endl;
		Ris2<<"------------------ END OF FILE ------------------"<<endl;
		Ris2.close();
}


// FORZA ELETTRICA - con moto associato
void elettr(double A, double B, double d){							// IN FUTURO CHIEDERE LA COSTANTE DIELETTRICA (NON SOLO EPSI)
	double Kappa=1/(4*pi*epsi);
	double vA=0, vB=0, pA=A*vA, pB=B*vB, F=(G*A*B)/pow(d,2.0), h=0;

	ofstream Ris3;
	Ris3.open("ForzaElettrica.txt", ios::out);

	Ris3<<"T [s]\td [m]\tp(A) [kg*m/s]\tp(B) [kg*m/s]\tF [N]\tvA [m/s]\tvB [m/s]"<<endl;

	for(double t=0; d>=0; t+=fq){
		Ris3<<t<<"\t"<<d<<"\t"<<pA<<"\t"<<pB<<"\t"<<F<<"\t"<<vA<<"\t"<<vB<<endl; // Stampa su file
		F=((Kappa*A*B)/pow(d,2.0));
		vA=(F/A)*t;   pA=A*vA;
		vB=(F/B)*t;   pB=B*vB;

		// ATTENZIONE !!! A e B NON SONO LE MASSE DEGLI OGGETTI, BENSI' LE LORO CARICHE (DA RICHIEDERE LA MASSA)
		h=(Kappa*A*B)/((0.5*A*pow(vA,2.0))+(0.5*B*pow(vB,2.0))-((Kappa*A*B)/d))+d;
		d+=h;
	}
	// Da inserire in futuro il richiamo alla procedura degli urti: quando i due corpi sono a distanza d=0 (oppure d=r1+r2)
	// Da richiedere prima se l'urto sara' elastico o meno, poi inviare alla procedura corrispondente
	cout<<"File salvato con successo."<<endl;
	Ris3<<"------------------ END OF FILE ------------------"<<endl;
	Ris3.close();
}


// URTI
void urtoE(double mA, double mB, double vA, double vB){
	double pA, pB, pAp, pBp, Ki, Kf, vBp=0, vAp=0;
	cout<<"NON IN FUNZIONE."<<endl;
	// p_i = p_f inoltre K_i = K_f
	// Conosco mA, mB e vA, vB

	pA = mA*vA;    pB = mB*vB;		 // Quantita' di moto iniziali
	pAp = mA*vAp;  pBp = mB*vBp;	 // Quantita' di moto finali

	Ki = 0.5*(mA*pow(vA, 2.0)+mB*pow(vB, 2.0));
	Kf = 0.5*(mA*pow(vAp, 2.0)+mB*pow(vBp, 2.0));
	
	cout<<"Simulazione terminata: controllare il file generato.";
	system("PAUSE");
}
void urtoA(double mA, double mB, double vA, double vB){
	cout<<"NON IN FUNZIONE."<<endl;
}


// PIANO INCLINATO
void PianoInclinato(double angle, double vi, double h){
	cout<<"NON IN FUNZIONE."<<endl;
}


// PRINCIPIO DI ARCHIMEDE
void Archy(){
	cout<<"NON IN FUNZIONE."<<endl; // da aggiungere se l'oggetto in caduta libera cade in acqua
}

// -------------------------------------------------------------

// PROCEDURA PER OTTENERE I DATI DA ELABORARE
void getDat(int ka){
	double x=0, y=0, teta=0, v=0, vB=0;
	int ko;
	cout<<"Per iniziare, inserire i seguenti dati.\n"<<endl;
	
	switch(ka){ // richiesta inserimento dei vari dati a seconda della simulazione scelta
		case 1:
			// - Richiedo inserimento dati per MOTO PARABOLICO
			cout<<"Altezza iniziale - slm [m]: "; cin>>y;
			cout<<"Angolo della velocita' rispetto all'orizzionte [gradi]: "; cin>>teta;
			cout<<"Modulo velocita' iniziale [m/s]: "; cin>>v;
			moti(x, y, v, teta*pi/180.0);
			break;

		case 2:
			// - Richiedo inserimento dati per URTI
			cout<<"Che tipo di urto? Elastico [1] o Anelastico [2]: ";			// NON VIENE CONSIDERATO L'ANGOLO DELLA VELOCITA'
			cin>>ko;															// DEI DUE CORPI !!!
			switch(ko){
				case 1:
					// Urto Elastico
					cout<<"Massa oggetto A [kg]: "; cin>>x;
					cout<<"Velocita' oggetto A [m/s]: "; cin>>v;
					cout<<"Massa oggetto B [kg]: "; cin>>y;
					cout<<"Velocita' oggetto B [m/s]: "; cin>>vB;
					urtoE(x, y, v, vB);
					break;
				case 2:
					// Urto Anelastico
					cout<<"Massa oggetto A [kg]: "; cin>>x;
					cout<<"Velocita' oggetto A [m/s]: "; cin>>v;
					cout<<"Massa oggetto B [kg]: "; cin>>y;
					cout<<"Velocita' oggetto B [m/s]: "; cin>>vB;
					urtoA(x, y, v, vB);
					break;
				default:
					cout<<"Errore!"<<endl;
					break;
			}
			break;

		case 3:
			// - Richiedo inserimento dati per FORZA ELETTRICA
			cout<<"Carica A [microC]: "; cin>>x;
			cout<<"Carica B [microC]: "; cin>>y;
			cout<<"Distanza fra A e B [m]: "; cin>>teta;
			x*=1e-6;
			y*=1e-6;
			elettr(x, y, teta);
			break;

		case 4:
			// - Richiedo inserimento dati per FORZA GRAVITAZIONALE
			cout<<"Massa oggetto A [kg]: "; cin>>x;
			cout<<"Massa oggetto B [kg]: "; cin>>y;
			cout<<"Distanza fra A e B [m]: "; cin>>teta;
			grav(x, y, teta);
			break;
	}
}

// PROCEDURA PER PROSEGUIRE UNA SIMULAZIONE IN UN ALTRO AMBIENTE (es. moto parabolico + rimbalzo pallina)
void prosDat(int scelta){
	cout<<"NON IN FUNZIONE."<<endl;
}

// -------------------------------------------------------------


int main(){ // menu del programma
	int ope;

	cout<<"Benvenuto nel simulatore. [v."<<versione<<"]"<<endl;
	cout<<"Programmato da Carlo Zambaldo, 5CS, Liceo G.Fracastoro, Verona IT"<<endl;
	cout<<"Per informazioni: carlo.zambaldo@gmail.com"<<endl<<endl;
	system("PAUSE");
	do{
		system("CLS");
		cout<<"Simulazioni possibili:\n1 - Moti in 2D (parabolico)\n2 - Urti \n3 - Forza Elettrica \n4 - Forza  Gravitazionale ";
		cout<<"\n5 - Piano Inclinato \n6 - Moto Gravita' + Principio di Archimede \n7 - Forze apparenti \n0 - Esci"<<endl;
		cout<<"\nInserire il numero dell'operazione da svolgere: ";
		cin>>ope;
		if((ope>=0)&&(ope<5)){
			if(ope!=0)
				getDat(ope);
		} else {
			cout<<"Errore!"<<endl;
		}
	} while(ope!=0);

	cout<<"Programma terminato."<<endl;
	cout<<"ATTENZIONE! - Pregasi di spostare i File creati prima della prossima apertura del programma."<<endl;
	system("PAUSE");
	return 0;
}
