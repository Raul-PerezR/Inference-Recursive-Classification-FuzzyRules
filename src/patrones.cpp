#include <cmath>
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <list>
#include <algorithm>    // std::sort
#include "patrones.h"
#include "example_set.h"

using namespace std;

void PausaP(){
	char ch;
	cout << "Pulsa una tecla para continuar" << endl;
	cin >> ch;
}

void deleteLine2(){
	putchar(27);    //codigo ascii decimal de ESC
	putchar('[');
	putchar('0');    //cero para borrar desde el cursor hasta el final
	putchar('K');

	putchar(27);    //codigo ascii decimal de ESC
	putchar('[');
	putchar('1');   // Se mueve una linea arriba
	putchar('A');

	putchar(27);    //codigo ascii decimal de ESC
	putchar('[');
	putchar('0');    //cero para borrar desde el cursor hasta el final
	putchar('K');
}

//-----------------------------------------------------------------------------------------------------
Pattern::Pattern(const Pattern & x){
  diccionario = x.diccionario;
  n_ejemplos = x.n_ejemplos;
  n_clases = x.n_clases;
}

//-----------------------------------------------------------------------------------------------------
Pattern &Pattern::operator=(const Pattern &x){
  if (this != &x){
    diccionario = x.diccionario;
    n_ejemplos = x.n_ejemplos;
    n_clases = x.n_clases;
  }
  return *this;
}

//-----------------------------------------------------------------------------------------------------
//Mejor etiqueta dada una variable y un ejemplo
int Pattern::BetterEtiqueta(const VectorVar &V, int variable, const example_set & E, int ejemplo){

    int l=V.SizeDomain(variable), et=0;
    double max=0;

    for (int etiqueta=0; etiqueta<l && max<1; etiqueta++){
        double aux = V.Adaptacion(E.Data(ejemplo,variable),variable,etiqueta);
        if (aux>max){
            max = aux;
            et = etiqueta;
        }
    }
    return et;
}


//-----------------------------------------------------------------------------------------------------
void ProcesarResultados(TestResult & result){
  double total_cubiertos = 0;
  double total_no_cubiertos = 0;
  double total = 0;

  result.acierto_global = 0;
  result.acierto_sinNoCubiertos = 0;
  result.porcentaje_nuevos_patrones = 0;
  result.error_intrinseco = 0;

  for (int i=0; i<result.acc.size(); i++){
    result.acierto_global+= result.acc[i];
    total_no_cubiertos += result.no_cubiertos[i];
    total_cubiertos += result.cubiertos[i];
  }

  total = total_no_cubiertos + total_cubiertos;

  if (total_cubiertos>0){
    result.acierto_sinNoCubiertos = 100.0*result.acierto_global / total_cubiertos;
  }

  if (total>0){
    result.acierto_global = 100.0*result.acierto_global / total;
    result.porcentaje_nuevos_patrones = 100.0*total_no_cubiertos / total;
  }
  result.error_intrinseco = 100.0 - result.acierto_sinNoCubiertos;
}

//-----------------------------------------------------------------------------------------------------

void PintaResultadosTraining(const TestResult & result){
  cout << "                 Acierto Global: " << result.acierto_global << " ----------------------" <<endl;
  for (int i=0; i<result.acc.size(); i++){
    if (result.no_cubiertos[i]+result.cubiertos[i]==0){
      cout << "                        Clase " << i <<": "  << 0 << endl;
    }
    else {
      cout << "                        Clase " << i <<": "  << 100.0*result.acc[i] / (result.no_cubiertos[i]+result.cubiertos[i]) << endl;
    }
  }
  cout << "               Error Intrinseco: " << result.error_intrinseco << " ----------------------" << endl;
  for (int i=0; i<result.acc.size(); i++){
    if (result.cubiertos[i]==0){
      cout << "                        Clase " << i <<": "  << 0 << endl;
    }
    else {
      cout << "                        Clase " << i <<": "  << 100.0 - 100.0*result.acc[i] / result.cubiertos[i] << endl;
    }
  }
}



//-----------------------------------------------------------------------------------------------------

void PintaResultadosTest(const TestResult & result, bool conClases){
  cout << "                 Acierto Global: " << result.acierto_global << " ----------------------" <<endl;
  if (conClases){
    for (int i=0; i<result.acc.size(); i++){
      if (result.no_cubiertos[i]+result.cubiertos[i]==0){
        cout << "                        Clase " << i <<": "  << 0 << endl;
      }
      else {
        cout << "                        Clase " << i <<": "  << 100.0*result.acc[i] / (result.no_cubiertos[i]+result.cubiertos[i]) << endl;
      }
    }
  }
  cout << "Acierto Global (Solo cubiertos): " << result.acierto_sinNoCubiertos << " ----------------------" << endl;
  if (conClases){
    for (int i=0; i<result.acc.size(); i++){
      if (result.cubiertos[i]==0){
        cout << "                        Clase " << i <<": "  << 0 << endl;
      }
      else {
        cout << "                        Clase " << i <<": "  << 100.0*result.acc[i] / result.cubiertos[i] << endl;
      }
    }
  }
  cout << "               Error Intrinseco: " << result.error_intrinseco << " ----------------------" << endl;
  if (conClases){
    for (int i=0; i<result.acc.size(); i++){
      if (result.cubiertos[i]==0){
        cout << "                        Clase " << i <<": "  << 0 << endl;
      }
      else {
        cout << "                        Clase " << i <<": "  << 100.0 - 100.0*result.acc[i] / result.cubiertos[i] << endl;
      }
    }
  }
  cout << "     Porcentaje Nuevos Patrones: " << result.porcentaje_nuevos_patrones << " ----------------------" << endl;
  if (conClases){
    for (int i=0; i<result.acc.size(); i++){
      if (result.no_cubiertos[i]+result.cubiertos[i]==0){
        cout << "                        Clase " << i <<": "  << 0  << endl;
      }
      else {
        cout << "                        Clase " << i <<": "  << 100.0 * result.no_cubiertos[i] / (result.no_cubiertos[i]+result.cubiertos[i])  << endl;
      }
    }
  }
}

//-----------------------------------------------------------------------------------------------------
void CalcularMediaTestResult(vector<TestResult> lista, TestResult & result){
  result.cubiertos.clear();
  result.no_cubiertos.clear();
  result.acc.clear();
  result.error_intrinseco_porClase.clear();

  for (int i=0; i<lista[0].acc.size(); i++){
    result.cubiertos.push_back(0);
    result.no_cubiertos.push_back(0);
    result.acc.push_back(0);
    result.error_intrinseco_porClase.push_back(0);
  }
  for (int j=0; j<lista.size(); j++){
    for (int i=0; i<lista[j].acc.size(); i++){
      result.cubiertos[i] += lista[j].cubiertos[i];
      result.no_cubiertos[i] += lista[j].no_cubiertos[i];
      result.acc[i] += lista[j].acc[i];
      result.error_intrinseco_porClase[i] += lista[j].error_intrinseco_porClase[i];
    }
  }

  ProcesarResultados(result);
}


//-----------------------------------------------------------------------------------------------------
void CalcularAgregacionResult(vector<TestResult> lista, TestResult & result){
  result.cubiertos.clear();
  result.no_cubiertos.clear();
  result.acc.clear();
  result.error_intrinseco_porClase.clear();

  for (int i=0; i<lista[0].acc.size(); i++){
    result.cubiertos.push_back(0);
    result.no_cubiertos.push_back(0);
    result.acc.push_back(0);
    result.error_intrinseco_porClase.push_back(0);
  }
  for (int j=0; j<lista.size(); j++){
    for (int i=0; i<lista[j].acc.size(); i++){
      result.cubiertos[i] += lista[j].cubiertos[i];
      result.no_cubiertos[i] += lista[j].no_cubiertos[i];
      result.acc[i] += lista[j].acc[i];
      result.error_intrinseco_porClase[i] += lista[j].error_intrinseco_porClase[i];
    }
  }

  ProcesarResultados(result);
}


//-----------------------------------------------------------------------------------------------------
void Pattern::ExtraerPatronesBasicos(const example_set &Es, const VectorVar &V, TestResult & result){

  // Crear la variable result
  int n_clases = V.SizeDomain(V.Consecuente());
  result.cubiertos.clear();
  result.no_cubiertos.clear();
  result.acc.clear();
  result.error_intrinseco_porClase.clear();

  for (int i=0; i<n_clases; i++){
    result.cubiertos.push_back(0);
    result.no_cubiertos.push_back(0);
    result.acc.push_back(0);
    result.error_intrinseco_porClase.push_back(0);
  }

  // Construyo el patron
  for (int i=0; i<Es.N_Examples(); i++){
    string aux="";
    for (int j=0; j<V.N_Antecedente(); j++){
      string auxaux="";
      for (int l=0; l<V.SizeDomain(j); l++){
        auxaux.push_back('0');
      }
      //auxaux.push_back(' ');
      auxaux[BetterEtiqueta(V,j,Es,i)]='1';
      aux += auxaux;
    }

    auto it = diccionario.find(aux);
    if (it==diccionario.end()){
      info d;
      for (int k=0; k<V.SizeDomain(V.Consecuente())+1; k++){
        d.conseq.push_back(0);
			}
      d.conseq[BetterEtiqueta(V,V.Consecuente(),Es,i)]++;
      d.eje_list[BetterEtiqueta(V,V.Consecuente(),Es,i)]++;
      diccionario.insert(pair<string,info > (aux,d));
    }
    else {
      it->second.conseq[BetterEtiqueta(V,V.Consecuente(),Es,i)]++;
			(it->second.eje_list[BetterEtiqueta(V,V.Consecuente(),Es,i)])++;
    }
  }



  auto it = diccionario.begin();
  double sumatotal=0, aciertostotal=0;
  while (it!=diccionario.end()){
    double suma=0;
    int mayor=0;
    for (int i=0; i<it->second.conseq.size(); i++){
      suma += it->second.conseq[i];
      result.cubiertos[i] += it->second.conseq[i];
      if (it->second.conseq[mayor]<it->second.conseq[i]){
        mayor = i;
      }
    }
    aciertostotal += it->second.conseq[mayor];
    result.acc[mayor] += it->second.conseq[mayor];
    sumatotal += suma;
    it++;
  }

  //cout << "Ejemplos: " << sumatotal << endl;
  //cout << "Aciertos: " << 100.0 * (aciertostotal / sumatotal) << endl;
  ProcesarResultados(result);

}



//-----------------------------------------------------------------------------------------------------
void Pattern::ExtraerPatronesBasicosOriginalWM(const example_set &Es, const VectorVar &V, TestResult & result){

  // Crear la variable result
  int n_clases = V.SizeDomain(V.Consecuente());
  result.cubiertos.clear();
  result.no_cubiertos.clear();
  result.acc.clear();
  result.error_intrinseco_porClase.clear();

  for (int i=0; i<n_clases; i++){
    result.cubiertos.push_back(0);
    result.no_cubiertos.push_back(0);
    result.acc.push_back(0);
    result.error_intrinseco_porClase.push_back(0);
  }

  // Construyo el patron
	double matching;
  for (int i=0; i<Es.N_Examples(); i++){
    string aux="";
		matching = 1.0;
    for (int j=0; j<V.N_Antecedente(); j++){
      string auxaux="";
      for (int l=0; l<V.SizeDomain(j); l++){
        auxaux.push_back('0');
      }
      //auxaux.push_back(' ');
			int localLabel = BetterEtiqueta(V,j,Es,i);
      auxaux[localLabel]='1';
			matching = matching * V.Adaptacion(Es.Data(i,j),j,localLabel);
      aux += auxaux;
    }

    auto it = diccionario.find(aux);
    if (it==diccionario.end()){
      info d;
			for (int k=0; k<V.SizeDomain(V.Consecuente())+1; k++){
        d.conseq.push_back(0);
				d.eje_list.push_back(0);
			}
      d.conseq[BetterEtiqueta(V,V.Consecuente(),Es,i)]+= matching;
      d.eje_list[BetterEtiqueta(V,V.Consecuente(),Es,i)]++;
      diccionario.insert(pair<string,info > (aux,d));
    }
    else {
      it->second.conseq[BetterEtiqueta(V,V.Consecuente(),Es,i)]+= matching;
			(it->second.eje_list[BetterEtiqueta(V,V.Consecuente(),Es,i)])++;
    }
  }



  auto it = diccionario.begin();
  double sumatotal=0, aciertostotal=0;
  while (it!=diccionario.end()){
    double suma=0;
    int mayor=0;
    for (int i=0; i < it->second.conseq.size(); i++){
      suma += it->second.eje_list[i];
      result.cubiertos[i] += it->second.eje_list[i];
      if (it->second.conseq[mayor] < it->second.conseq[i]){
        mayor = i;
      }
    }
    aciertostotal += it->second.eje_list[mayor];
    result.acc[mayor] += it->second.eje_list[mayor];
    sumatotal += suma;
    it++;
  }

	cout << "Patrones: " << diccionario.size() << endl;
	/*char ch;
	cout << "pulsa una tecla en extraer patrones originales " << endl;
	cin >> ch;*/

  //cout << "Ejemplos: " << sumatotal << endl;
  //cout << "Aciertos: " << 100.0 * (aciertostotal / sumatotal) << endl;
  ProcesarResultados(result);

}

void ponerEnCadena(string &cadena, int pos, string subcadena){
	for (int i=0; i<subcadena.size(); i++)
	  cadena[pos+i] =subcadena[i];
}

bool Ordenar_infoUp (const infoUp &x, const infoUp &y){
	return (x.key > y.key);
}


//-----------------------------------------------------------------------------------------------------
void Pattern::Aprendizaje_RecursivoUnEjemplo_WM_TFM_Ruben(const example_set &Es, const VectorVar &V, const int eje,
		 string cadena, int actualVar, double adapt, const vector<infoUp> &trozos, int &current_tries, int number_tries){
  int j = actualVar;

	if (adapt<=0){  // la adaptacion con el ejemplo es 0
		return;
	}
  else if (current_tries >= number_tries){ // Condicion de parada indica mala regla
    return;
  }
  else if (j == V.N_Antecedente()){ // Condicion de parada cuando se ha completado una regla
		current_tries++;
		//cout << current_tries << ") " << cadena << " -->" << adapt;
		//cout << endl;
		auto it = diccionario.find(cadena);
		info d;
    if (it==diccionario.end()){
			for (int k=0; k<V.SizeDomain(V.Consecuente())+1; k++){
        d.conseq.push_back(0);
				d.eje_list.push_back(0);
			}
      d.conseq[BetterEtiqueta(V,V.Consecuente(),Es,eje)]+= adapt;
      d.eje_list[BetterEtiqueta(V,V.Consecuente(),Es,eje)]++;
      diccionario.insert(pair<string,info > (cadena,d));
			/*cout << "Tama Dicionario: " << diccionario.size() << endl;
			cout << "Se añade " << cadena << " [";
			for (int k=0; k<V.SizeDomain(V.Consecuente())-1; k++){
				cout << d.conseq[k] << " |";
			}
			cout << d.conseq[V.SizeDomain(V.Consecuente())] << "]\n";*/

    }
    else {
      it->second.conseq[BetterEtiqueta(V,V.Consecuente(),Es,eje)]+= adapt;
			(it->second.eje_list[BetterEtiqueta(V,V.Consecuente(),Es,eje)])++;
			/*cout << "-------> " << cadena << " [";
			for (int k=0; k<V.SizeDomain(V.Consecuente())-1; k++){
				cout << it->second.conseq[k] << " |";
			}
			cout << it->second.conseq[V.SizeDomain(V.Consecuente())] << "]\n";
			char ch;
			cin >> ch;*/
    }
  }
  else { // Cuerpo de la recursion
    for (int i=0; i<trozos[j].trozoAcambiar.size();i++){
			string copy_cadena = cadena;
			ponerEnCadena(copy_cadena,trozos[j].posCadenainicio,trozos[j].trozoAcambiar[i].subCadena);
			// T-Norm product
			//cout << copy_cadena << endl;
			Aprendizaje_RecursivoUnEjemplo_WM_TFM_Ruben(Es,V,eje,copy_cadena, j+1, adapt * trozos[j].trozoAcambiar[i].adaptacion, trozos, current_tries, number_tries);
		}
  }
}


//-----------------------------------------------------------------------------------------------------
void Pattern::ExtraerPatronesBasicosAproximacionTFMRuben(const example_set &E, const VectorVar &V, TestResult & result, int tries){

  cout << "Estoy en ExtraerPatronesBasicosAproximacionTFMRuben\n" << endl;
  // Crear la variable result
  int n_clases = V.SizeDomain(V.Consecuente());
  result.cubiertos.clear();
  result.no_cubiertos.clear();
  result.acc.clear();
  result.error_intrinseco_porClase.clear();

  for (int i=0; i<n_clases; i++){
    result.cubiertos.push_back(0);
    result.no_cubiertos.push_back(0);
    result.acc.push_back(0);
    result.error_intrinseco_porClase.push_back(0);
  }

	vector<infoUp>trozos;
  // Construyo el patron
  for (int i=0; i<E.N_Examples(); i++){
		//pinto el ejemplos
		/*cout << "(" << E.Data(i,0);
		for (int j=1; j<V.N_Antecedente(); j++){
			cout <<", " << E.Data(i,j);
		}
		cout << ") ----- >" <<  E.Data(i,V.Consecuente()) << endl;*/
		deleteLine2();
		cout << "Learning:   " << (i*100.0)/E.N_Examples() << "\%"<<endl;
		int current_tries = 0;
		string coded_rule="";
    trozos.clear();

		int tama=0;
		for (int j=0; j<V.N_Antecedente(); j++){
			string coded_var="";
			infoUp datos;

			datos.variable = j;
			datos.posCadenainicio = tama;
			datos.key = 1;

			for (int l=0; l<V.SizeDomain(j); l++){
				coded_rule.push_back('-');
				coded_var.push_back('0');
			}

			for (int l=0; l<V.SizeDomain(j); l++){
				info2 datos2;
				double d = V.Adaptacion(E.Data(i,j),j,l);
				if (d>0){
					datos2.adaptacion = d;
					coded_var[l]='1';
					datos2.subCadena = coded_var;
					coded_var[l]='0';
					// meter las modificaciones ordenadas.
					if (datos.trozoAcambiar.size() != 1)
						datos.trozoAcambiar.push_back(datos2);
					else {
						if (datos.trozoAcambiar[0].adaptacion > datos2.adaptacion){
							datos.key = datos.trozoAcambiar[0].adaptacion - datos2.adaptacion;
							datos.trozoAcambiar.push_back(datos2);
						}
						else {
							datos.key = datos2.adaptacion - datos.trozoAcambiar[0].adaptacion;
							datos.trozoAcambiar.push_back(datos2);
							swap(datos.trozoAcambiar[0],datos.trozoAcambiar[1]);
						}
					}
				}
			}
      // Aniado eliminar esa variables
			/*info2 datos3;
			datos3.adaptacion = 1.0;
			for (int l=0; l<V.SizeDomain(j); l++)
				coded_var[l]='1';
			datos3.subCadena = coded_var;
			datos.trozoAcambiar.push_back(datos3);*/


			trozos.push_back(datos);
			tama+=V.SizeDomain(j);
		}


	  //Incluir el método sort (requiere #include <algorithm>) para ordenar el vector.
		/*sort(trozos.begin(), trozos.end(), Ordenar_infoUp);
		for (int i =0; i < trozos.size(); i++){
			cout << "variable: " << trozos[i].variable << "\t" << trozos[i].key;
			for (int j = 0; j< trozos[i].trozoAcambiar.size(); j++){
				cout << "   (" << trozos[i].trozoAcambiar[j].subCadena << ")"  << trozos[i].trozoAcambiar[j].adaptacion;
			}
			cout << endl;
		}
	  cout << endl;*/


		Aprendizaje_RecursivoUnEjemplo_WM_TFM_Ruben(E,V,i,coded_rule, 0, 1, trozos, current_tries, tries);
		/*char ch;
		cout << "pulsa una tecla para continuar \n";
		cin >> ch;*/
  }




  // Resumen final de resultados
  auto it = diccionario.begin();
  double sumatotal=0, aciertostotal=0;
  while (it!=diccionario.end()){
    double suma=0, suma_adapt=0;
    int mayor=0;
    for (int i=0; i < it->second.conseq.size()-1; i++){
      suma += it->second.eje_list[i];
			suma_adapt += it->second.conseq[i];
      result.cubiertos[i] += it->second.eje_list[i];
      if (it->second.conseq[mayor] < it->second.conseq[i]){
        mayor = i;
      }
    }
    aciertostotal += it->second.eje_list[mayor];
    result.acc[mayor] += it->second.eje_list[mayor];
    sumatotal += suma;
		it->second.eje_list[it->second.conseq.size()-1] = suma;
    /******** salida de los antecedentes encontrados **********/
		/*cout << it->first << "\t(" << it->second.eje_list[0];
		for (int i=1; i < it->second.conseq.size(); i++){
			cout << " | " << it->second.eje_list[i];
		}
		cout << ")";
		cout << endl;*/
		/**********************************************************/
		/******** salida de los antecedentes porcentaje acumulacion adaptacion **********/
		/*int a = (it->second.conseq[0] / suma_adapt)*100;
		cout << it->first << "\t(" << a;
		for (int i=1; i < it->second.conseq.size(); i++){
			a = (it->second.conseq[i] / suma_adapt)*100;
			cout << " | " << a;
		}
		cout << ")";
		cout << endl;*/
		/**********************************************************/

    it++;
  }

	cout << "Patrones: " << diccionario.size() << endl;
	/*char ch;
	cout << "pulsa una tecla" << endl;
	cin >> ch;*/

  //cout << "Ejemplos: " << sumatotal << endl;
  //cout << "Aciertos: " << 100.0 * (aciertostotal / sumatotal) << endl;
  ProcesarResultados(result);

}


//-----------------------------------------------------------------------------------------------------

double min(double a, double b){
  if (a<=b) return a; else return b;
}





void OrdenarElVector (vector<infoUp> & v){
	for (int i=0; i< v.size()-1; i++){
		for (int j=i+1; j<v.size(); j++){
			if (v[i].key>v[j].key){
				infoUp aux = v[i];
				v[i] = v[j];
				v[j] = aux;
			}
		}
	}
}

//-----------------------------------------------------------------------------------------------------

int Pattern::InferenciaRecursivaOptimalizada(const example_set &E, const VectorVar &V, const int eje, bool PCF, double &mu, string &antecedenteSeleccionado, int TriesNumber){
	string coded_rule="";
	vector<infoUp>trozos;


	int tama=0;
	for (int j=0; j<V.N_Antecedente(); j++){
		string coded_var="";
		infoUp datos;

		datos.variable = j;
		datos.posCadenainicio = tama;
		datos.key = 1;

		for (int l=0; l<V.SizeDomain(j); l++){
			coded_rule.push_back('0');
			coded_var.push_back('0');
		}

		for (int l=0; l<V.SizeDomain(j); l++){
			info2 datos2;
			double d = V.Adaptacion(E.Data(eje,j),j,l);
			if (d>0){
				datos2.adaptacion = d;
				coded_var[l]='1';
				datos2.subCadena = coded_var;
				coded_var[l]='0';
				// meter las modificaciones ordenadas.
				if (datos.trozoAcambiar.size() != 1)
					datos.trozoAcambiar.push_back(datos2);
				else {
					if (datos.trozoAcambiar[0].adaptacion > datos2.adaptacion){
						datos.key = datos.trozoAcambiar[0].adaptacion - datos2.adaptacion;
						datos.trozoAcambiar.push_back(datos2);
					}
					else {
						datos.key = datos2.adaptacion - datos.trozoAcambiar[0].adaptacion;
						datos.trozoAcambiar.push_back(datos2);
						swap(datos.trozoAcambiar[0],datos.trozoAcambiar[1]);
					}
				}
			}
		}
		trozos.push_back(datos);
		tama+=V.SizeDomain(j);
	}


  //Incluir el método sort (requiere #include <algorithm>) para ordenar el vector.
	sort(trozos.begin(), trozos.end(), Ordenar_infoUp);

	/*for (int i =0; i < trozos.size(); i++){
		cout << "variable: " << trozos[i].variable << "\t" << trozos[i].key;
		for (int j = 0; j< trozos[i].trozoAcambiar.size(); j++){
			cout << "   (" << trozos[i].trozoAcambiar[j].subCadena << ")"  << trozos[i].trozoAcambiar[j].adaptacion;
		}
		cout << endl;
	}
  cout << endl;*/

	int clase=-1;
	double umbral=0.0;
	TestearRecursivoUnEjemploEficente(E,V,eje,coded_rule,0,1.0,umbral,clase,PCF,trozos, antecedenteSeleccionado, TriesNumber);
	mu=umbral;
	return clase;
}








//-----------------------------------------------------------------------------------------------------
void Pattern::TestearRecursivoUnEjemploEficente(const example_set &E, const VectorVar &V, const int eje,
		 string cadena, int actualVar, double adapt, double &umbral, int &clase, bool PCF, const vector<infoUp> &trozos, string &antecedenteSeleccionado, int &TriesNumber){
  int j = actualVar;
  //deleteLine2();
	//if (j>=0 and j<V.N_Antecedente()){
  //	cout << trozos[j].variable << "\t" << cadena << " A:" << adapt << " U:" << umbral << endl;
	//}
	//else {
	//	cout << "\t" << cadena << " A:" << adapt << " U:" << umbral << endl;
	//}
	//char ch;
	//cin >> ch;

  if (TriesNumber<=0) // Condicion de parada por superar el número máximo de reglas a evaluar
	  return;
  else if (adapt <= umbral){ // Condicion de parada indica mala regla
				if (clase != -1) TriesNumber--;
    return;
  }
  else if (j == V.N_Antecedente()){ // Condicion de parada cuando se ha completado una regla
		TriesNumber--;
		//cout << "Intentos: " << TriesNumber << "  Umbral: " << umbral << " Clase: " << clase << endl << endl;
    auto it = diccionario.find(cadena);
    int mayor;
		double ejemplos;
    if (it!=diccionario.end()){
      ejemplos = it->second.conseq[0];
      mayor=0;
      for (int k=1; k<V.SizeDomain(V.Consecuente()); k++){
        ejemplos = ejemplos + it->second.conseq[k];
        if (it->second.conseq[k]>it->second.conseq[mayor])
        	mayor=k;
      }
      double peso;

			if (PCF){
				peso = (1.0 * it->second.conseq[mayor] - (ejemplos-it->second.conseq[mayor])) / (ejemplos);
			}
			else {
				peso = (1.0 * it->second.conseq[mayor])/ (ejemplos);
			}
      adapt = adapt * peso;
      if (adapt>umbral){
        umbral = adapt;
        clase = mayor;
				antecedenteSeleccionado = cadena;
        //deleteLine2();
        //cout << "---->\t" << cadena << " T: " << umbral << " Cl: " << clase << " ClReal: " << E.Data(eje,V.Consecuente())<<endl;
				//cout << "---->\t" << cadena << " T: " << umbral << " Cl: " << clase << " ClReal: " << E.Data(eje,V.Consecuente())<< " peso: " << peso << endl;
				//char ch;
				//cout << "Pulsa Una Tecla ...." << endl;
				//cin >> ch;
      }
    }

  }
  else { // Cuerpo de la recursion
    for (int i=0; i<trozos[j].trozoAcambiar.size();i++){
			string copy_cadena = cadena;
			ponerEnCadena(copy_cadena,trozos[j].posCadenainicio,trozos[j].trozoAcambiar[i].subCadena);
			// T-Norm minimo
			//TestearRecursivoUnEjemploEficente(E,V,eje,copy_cadena, j+1, min(adapt, trozos[j].trozoAcambiar[i].adaptacion), umbral, clase, PCF, trozos);
			// T-Norm product
			TestearRecursivoUnEjemploEficente(E,V,eje,copy_cadena, j+1, adapt * trozos[j].trozoAcambiar[i].adaptacion, umbral, clase, PCF, trozos, antecedenteSeleccionado, TriesNumber);
			//deleteLine2();
		}
  }
}



//-----------------------------------------------------------------------------------------------------

int Pattern::InferenciaRecursivaOptimalizadaConProfundidadLimitada(const example_set &E, const VectorVar &V, const int eje, bool PCF, double &mu, string &antecedenteSeleccionado, int prof){
	string coded_rule="";
	vector<infoUp>trozos;


	int tama=0;
	for (int j=0; j<V.N_Antecedente(); j++){
		string coded_var="";
		infoUp datos;

		datos.variable = j;
		datos.posCadenainicio = tama;
		datos.key = 1;

		for (int l=0; l<V.SizeDomain(j); l++){
			coded_rule.push_back('0');
			coded_var.push_back('0');
		}

		for (int l=0; l<V.SizeDomain(j); l++){
			info2 datos2;
			double d = V.Adaptacion(E.Data(eje,j),j,l);
			if (d>0){
				datos2.adaptacion = d;
				coded_var[l]='1';
				datos2.subCadena = coded_var;
				coded_var[l]='0';
				// meter las modificaciones ordenadas.
				if (datos.trozoAcambiar.size() != 1)
					datos.trozoAcambiar.push_back(datos2);
				else {
					if (datos.trozoAcambiar[0].adaptacion > datos2.adaptacion){
						datos.key = datos.trozoAcambiar[0].adaptacion - datos2.adaptacion;
						datos.trozoAcambiar.push_back(datos2);
					}
					else {
						datos.key = datos2.adaptacion - datos.trozoAcambiar[0].adaptacion;
						datos.trozoAcambiar.push_back(datos2);
						swap(datos.trozoAcambiar[0],datos.trozoAcambiar[1]);
					}
				}
			}
		}
		trozos.push_back(datos);
		tama+=V.SizeDomain(j);
	}


  //Incluir el método sort (requiere #include <algorithm>) para ordenar el vector.
	sort(trozos.begin(), trozos.end(), Ordenar_infoUp);

	/*for (int i =0; i < trozos.size(); i++){
		cout << "variable: " << trozos[i].variable << "\t" << trozos[i].key;
		for (int j = 0; j< trozos[i].trozoAcambiar.size(); j++){
			cout << "   (" << trozos[i].trozoAcambiar[j].subCadena << ")"  << trozos[i].trozoAcambiar[j].adaptacion;
		}
		cout << endl;
	}
  cout << endl;*/

	int clase=-1;
	double umbral = 0.0;
	TestearRecursivoUnEjemploEficenteConProfundidadLimitada(E,V,eje,coded_rule,0,1.0,umbral,clase,PCF,trozos, antecedenteSeleccionado, prof);
	mu=umbral;
	return clase;
}


//-----------------------------------------------------------------------------------------------------
void Pattern::TestearRecursivoUnEjemploEficenteConProfundidadLimitada(const example_set &E, const VectorVar &V, const int eje,
		 string cadena, int actualVar, double adapt, double &umbral, int &clase, bool PCF, const vector<infoUp> &trozos, string &antecedenteSeleccionado, int prof){
  int j = actualVar;
  //deleteLine2();
	//if (j>=0 and j<V.N_Antecedente()){
  //	cout << trozos[j].variable << "\t" << cadena << " A:" << adapt << " U:" << umbral << endl;
	//}
	//else {
	//	cout << "\t" << cadena << " A:" << adapt << " U:" << umbral << endl;
	//}
	//char ch;
	//cin >> ch;

  if (prof < 0){ //Se alcanzó la profundidad maxima
		return;
	}
  else if (adapt <= umbral){ // Condicion de parada indica mala regla
    return;
  }
  else if (j == V.N_Antecedente()){ // Condicion de parada cuando se ha completado una regla
    auto it = diccionario.find(cadena);
    int mayor;
		double ejemplos;
    if (it!=diccionario.end()){
      ejemplos = it->second.conseq[0];
      mayor=0;
      for (int k=1; k<V.SizeDomain(V.Consecuente()); k++){
        ejemplos = ejemplos + it->second.conseq[k];
        if (it->second.conseq[k]>it->second.conseq[mayor])
        	mayor=k;
      }
      double peso;

			if (PCF){
				peso = (1.0 * it->second.conseq[mayor] - (ejemplos-it->second.conseq[mayor])) / (ejemplos);
			}
			else {
				peso = (1.0 * it->second.conseq[mayor])/ (ejemplos);
			}
      adapt = adapt * peso;
      if (adapt>umbral){
        umbral = adapt;
        clase = mayor;
				antecedenteSeleccionado = cadena;
        //deleteLine2();
        //cout << "---->\t" << cadena << " T: " << umbral << " Cl: " << clase << " ClReal: " << E.Data(eje,V.Consecuente())<<endl;
				//cout << "---->\t" << cadena << " T: " << umbral << " Cl: " << clase << " ClReal: " << E.Data(eje,V.Consecuente())<< " peso: " << peso << endl;
				//char ch;
				//cout << "Pulsa Una Tecla ...." << endl;
				//cin >> ch;
      }
    }

  }
  else { // Cuerpo de la recursion
    for (int i=0; i<trozos[j].trozoAcambiar.size();i++){
			string copy_cadena = cadena;
			ponerEnCadena(copy_cadena,trozos[j].posCadenainicio,trozos[j].trozoAcambiar[i].subCadena);
			// T-Norm minimo
			//TestearRecursivoUnEjemploEficente(E,V,eje,copy_cadena, j+1, min(adapt, trozos[j].trozoAcambiar[i].adaptacion), umbral, clase, PCF, trozos);
			// T-Norm product
			TestearRecursivoUnEjemploEficenteConProfundidadLimitada(E,V,eje,copy_cadena, j+1, adapt * trozos[j].trozoAcambiar[i].adaptacion, umbral, clase, PCF, trozos, antecedenteSeleccionado, prof-i);
			//deleteLine2();
		}
  }
}




//-----------------------------------------------------------------------------------------------------

int Pattern::InferenciaRecursivaOptimalizadaConProfundidadLimitada_SalidaPorProfundidades(const example_set &E, const VectorVar &V, const int eje, bool PCF, double &mu, string &antecedenteSeleccionado, int prof){
	string coded_rule="";
	vector<infoUp>trozos;


	int tama=0;
	for (int j=0; j<V.N_Antecedente(); j++){
		string coded_var="";
		infoUp datos;

		datos.variable = j;
		datos.posCadenainicio = tama;
		datos.key = 1;

		for (int l=0; l<V.SizeDomain(j); l++){
			coded_rule.push_back('0');
			coded_var.push_back('0');
		}

		for (int l=0; l<V.SizeDomain(j); l++){
			info2 datos2;
			double d = V.Adaptacion(E.Data(eje,j),j,l);
			if (d>0){
				datos2.adaptacion = d;
				coded_var[l]='1';
				datos2.subCadena = coded_var;
				coded_var[l]='0';
				// meter las modificaciones ordenadas.
				if (datos.trozoAcambiar.size() != 1)
					datos.trozoAcambiar.push_back(datos2);
				else {
					if (datos.trozoAcambiar[0].adaptacion > datos2.adaptacion){
						datos.key = datos.trozoAcambiar[0].adaptacion - datos2.adaptacion;
						datos.trozoAcambiar.push_back(datos2);
					}
					else {
						datos.key = datos2.adaptacion - datos.trozoAcambiar[0].adaptacion;
						datos.trozoAcambiar.push_back(datos2);
						swap(datos.trozoAcambiar[0],datos.trozoAcambiar[1]);
					}
				}
			}
		}
		trozos.push_back(datos);
		tama+=V.SizeDomain(j);
	}


  //Incluir el método sort (requiere #include <algorithm>) para ordenar el vector.
	sort(trozos.begin(), trozos.end(), Ordenar_infoUp);

	/*for (int i =0; i < trozos.size(); i++){
		cout << "variable: " << trozos[i].variable << "\t" << trozos[i].key;
		for (int j = 0; j< trozos[i].trozoAcambiar.size(); j++){
			cout << "   (" << trozos[i].trozoAcambiar[j].subCadena << ")"  << trozos[i].trozoAcambiar[j].adaptacion;
		}
		cout << endl;
	}
  cout << endl;*/

	vector<int> clase;
	vector<double> umbral;
	vector<pair<string,double> > estructura;
	pair<string,double> auxpair;
	auxpair.first = "-";
	auxpair.second =-1;
	for (int p=0; p<prof; p++){
		clase.push_back(-1);
		umbral.push_back(0.0);
		estructura.push_back(auxpair);
	}

	TestearRecursivoUnEjemploEficenteConProfundidadLimitada_SalidaPorProfundidades(E,V,eje,coded_rule,0,1.0,umbral,clase,PCF,trozos, antecedenteSeleccionado, prof, estructura);

  int decidirClase = -1;
	for (int p=0; p<prof; p++){
		if (decidirClase == -1) decidirClase = clase[p];
		cout << "umbral " << p << ": " << umbral[p] << " C: " << clase[p] << " -->\t" << estructura[p].first << " , " << estructura[p].second << endl;
	}
	PausaP();



	return decidirClase;
}



//-----------------------------------------------------------------------------------------------------
void Pattern::TestearRecursivoUnEjemploEficenteConProfundidadLimitada_SalidaPorProfundidades(const example_set &E, const VectorVar &V, const int eje,
		 string cadena, int actualVar, double adapt, vector<double> &umbral, vector<int> &clase, bool PCF, const vector<infoUp> &trozos, string &antecedenteSeleccionado,
		 int prof, vector<pair<string,double> > &estructura){
  int j = actualVar;
  //deleteLine2();
	//if (j>=0 and j<V.N_Antecedente()){
  //	cout << trozos[j].variable << "\t" << cadena << " A:" << adapt << " U:" << umbral << endl;
	//}
	//else {
	//	cout << "\t" << cadena << " A:" << adapt << " U:" << umbral << endl;
	//}
	//char ch;
	//cin >> ch;

  if (prof < 0){ //Se alcanzó la profundidad maxima
		return;
	}
  else if (adapt <= umbral[prof]){ // Condicion de parada indica mala regla
    return;
  }
  else if (j == V.N_Antecedente()){ // Condicion de parada cuando se ha completado una regla
    auto it = diccionario.find(cadena);
    int mayor;
		double ejemplos;
    if (it!=diccionario.end()){
      ejemplos = it->second.conseq[0];
      mayor=0;
      for (int k=1; k<V.SizeDomain(V.Consecuente()); k++){
        ejemplos = ejemplos + it->second.conseq[k];
        if (it->second.conseq[k]>it->second.conseq[mayor])
        	mayor=k;
      }
      double peso;

			if (PCF){
				peso = (1.0 * it->second.conseq[mayor] - (ejemplos-it->second.conseq[mayor])) / (ejemplos);
			}
			else {
				peso = (1.0 * it->second.conseq[mayor])/ (ejemplos);
			}

      adapt = adapt * peso;
      if (adapt>umbral[prof]){
        umbral[prof] = adapt;
        clase[prof] = mayor;
				antecedenteSeleccionado = cadena;
        //deleteLine2();
        //cout << "---->\t" << cadena << " T: " << umbral << " Cl: " << clase << " ClReal: " << E.Data(eje,V.Consecuente())<<endl;
				//cout << "---->\t" << cadena << " T: " << umbral << " Cl: " << clase << " ClReal: " << E.Data(eje,V.Consecuente())<< " peso: " << peso << endl;
				//char ch;
				//cout << "Pulsa Una Tecla ...." << endl;
				//cin >> ch;
				// Se aniade a la estructura de decisiones
				pair<string,double> unpar;
				unpar.first = cadena;
				unpar.second = (1.0 * it->second.conseq[mayor])/ (ejemplos); // El peso
				estructura[prof] = unpar;

      }
    }

  }
  else { // Cuerpo de la recursion
    for (int i=0; i<trozos[j].trozoAcambiar.size();i++){
			string copy_cadena = cadena;
			ponerEnCadena(copy_cadena,trozos[j].posCadenainicio,trozos[j].trozoAcambiar[i].subCadena);
			// T-Norm minimo
			//TestearRecursivoUnEjemploEficente(E,V,eje,copy_cadena, j+1, min(adapt, trozos[j].trozoAcambiar[i].adaptacion), umbral, clase, PCF, trozos);
			// T-Norm product
			TestearRecursivoUnEjemploEficenteConProfundidadLimitada_SalidaPorProfundidades(E,V,eje,copy_cadena, j+1, adapt * trozos[j].trozoAcambiar[i].adaptacion, umbral, clase, PCF, trozos, antecedenteSeleccionado, prof-i, estructura);
			//deleteLine2();
		}
  }
}




//-----------------------------------------------------------------------------------------------------
void Pattern::TestearRecursivoUnEjemplo(const example_set &E, const VectorVar &V, const int eje, string cadena, int actualVar, double adapt, double &umbral, int &clase, bool PCF, string &antecedenteSeleccionado, int &TriesNumber){
  int j = actualVar;
  //deleteLine2();
  //cout << j << "\t" << cadena << " A:" << adapt << " U:" << umbral << endl;

  if (TriesNumber<=0){  // Condicion de parada por haber evaluado el numero maximo de reglas
		return;
	}
  else if (adapt <= umbral){ // Condicion de parada indica mala regla
		if (clase != -1) TriesNumber--;
    return;
  }
  else if (j == V.N_Antecedente()){ // Condicion de parada cuando se ha completado una regla
		TriesNumber--;
    auto it = diccionario.find(cadena);
    int mayor;
    if (it!=diccionario.end()){
      double ejemplos = it->second.conseq[0];
      mayor=0;
      for (int k=1; k<n_clases; k++){
        ejemplos += it->second.conseq[k];
        if (it->second.conseq[k]>it->second.conseq[mayor])
        mayor=k;
      }
      double peso;

			if (PCF){
				peso = (1.0 * it->second.conseq[mayor] - (ejemplos-it->second.conseq[mayor])) / (ejemplos);
			}
			else {
				peso = (1.0 * it->second.conseq[mayor])/ (ejemplos);
			}
      adapt = adapt * peso;
      if (adapt>umbral){
        umbral = adapt;
        clase = mayor;
				antecedenteSeleccionado = cadena;
        //deleteLine2();
        //cout << "---->\t" << cadena << " T: " << umbral << " Cl: " << clase << " ClReal: " << E.Data(eje,V.Consecuente())<<endl;
				//char ch;
				//cout << "Pulsa Una Tecla ...." << endl;
				//cin >> ch;
      }
    }

  }
  else { // Cuerpo de la recursion
    string auxaux="";
    pair<double,string> pareja;
    std::map<double,string, std::greater<double> > lista;

    for (int l=0; l<V.SizeDomain(j); l++){
      auxaux.push_back('0');
    }

    for (int l=0; l<V.SizeDomain(j); l++){
      double d = V.Adaptacion(E.Data(eje,j),j,l);
      if (d>0){
        pareja.first = d;
        auxaux[l]='1';
        pareja.second = auxaux;
        auxaux[l]='0';
        lista.insert( pareja );
      }
    }

    for (auto it = lista.begin(); it != lista.end(); it++){
      //cout << endl;
			// T-Norm minimo
      //TestearRecursivoUnEjemplo(E,V,eje,cadena+it->second, j+1, min(adapt,it->first), umbral, clase);

			// T-Norm producto
			TestearRecursivoUnEjemplo(E,V,eje,cadena+it->second, j+1, adapt * it->first, umbral, clase, PCF, antecedenteSeleccionado, TriesNumber);
      //deleteLine2();
    }

  }
}

//-----------------------------------------------------------------------------------------------------
void Pattern::TestearRecursivo(const example_set &Es, const VectorVar &V, TestResult & result, bool PCF, int TriesNumber){

  int n_clases = V.SizeDomain(V.Consecuente());
  result.cubiertos.clear();
  result.no_cubiertos.clear();
  result.acc.clear();
  result.error_intrinseco_porClase.clear();

  for (int i=0; i<n_clases; i++){
    result.cubiertos.push_back(0);
    result.no_cubiertos.push_back(0);
    result.acc.push_back(0);
    result.error_intrinseco_porClase.push_back(0);
  }

  double aciertos = 0, noCub = 0;
  cout << endl;
  int old_value=0, new_value;
  for (int i=0; i<Es.N_Examples(); i++){
    string regla="";
    double umbral = 0;
    int clase_predicha = -1;
    new_value = 100*i/Es.N_Examples();
    if (new_value != old_value){
      deleteLine2();
      cout << "(%) Ejemplo " << 1.0*new_value << endl;
      old_value = new_value;
    }
		string antecedenteSeleccionado;
    TestearRecursivoUnEjemplo(Es,V,i,regla,0,1,umbral,clase_predicha, PCF, antecedenteSeleccionado, TriesNumber);
    //char ch;
    //cin >> ch;
    int clase_real = Es.Data(i,V.Consecuente());
    result.cubiertos[clase_real]++;
    if (clase_predicha==clase_real){
      result.acc[clase_real]++;
      aciertos++;
    }
    if (clase_predicha == -1){
      result.no_cubiertos[clase_real]++;
      noCub++;
    }
  }
  cout << "Aciertos: " << 100*aciertos/Es.N_Examples() << endl;
  cout << "Aciertos no cubiertos: " << 100*aciertos/(Es.N_Examples()-noCub) << endl;
  cout << "Porcent no cubiertos: " << 100*noCub/Es.N_Examples() << endl;


  ProcesarResultados(result);
}
//-----------------------------------------------------------------------------------------------------
int DistanciaHamming(string r1, string r2){
	int df=0;
	for (int i=0; i < r1.length(); i++)
	  if (r1[i]!=r2[i])
		  df++;
	return df;
}

//-----------------------------------------------------------------------------------------------------
int DistanciaHamming(string r1, string r2, int menosQueEsta){
	int df=0;
	for (int i=0; df <= menosQueEsta and i < r1.length(); i++)
	  if (r1[i]!=r2[i])
		  df++;
	return df;
}

//-----------------------------------------------------------------------------------------------------
int Pattern::TesteoDistanciaHamming(const example_set &Es, const VectorVar &V, const int i){
  // Construyo la regla que mejor se adapta al ejemplo i
	string r="";
	for (int j=0; j<V.N_Antecedente(); j++){
		string auxaux="";
		for (int l=0; l<V.SizeDomain(j); l++){
			auxaux.push_back('0');
		}
		//auxaux.push_back(' ');
		auxaux[BetterEtiqueta(V,j,Es,i)]='1';
		r += auxaux;
	}

  // Busco y evaluo las reglas con menor distancia hamming
	int distancia = r.length()+1;
	int reglasImplicadas = 0;
	vector<int> acumulado;
	for (int i = 0; i< n_clases; i++){
		acumulado.push_back(0);
	}

	for (auto it = diccionario.begin(); it != diccionario.end(); it++){
		int distTor1 = DistanciaHamming(r,it->first, distancia);
		if (distTor1 < distancia){
			reglasImplicadas = 1;
			distancia = distTor1;
			for (int i = 0; i< n_clases; i++){
				acumulado[i] = it->second.conseq[i];
			}
		}
		else if (distTor1 == distancia){
			reglasImplicadas++;
			for (int i = 0; i< n_clases; i++){
				acumulado[i] += it->second.conseq[i];
			}
		}
	}

	//cout << "Dist: " << distancia << "  ReglasImplicadas: " << reglasImplicadas << " ";
	int mejor_clase = -1, mejor_actual = 0;
	for (int i = 0; i< n_clases; i++){
		//cout << acumulado[i] << " ";
		if (acumulado[i]> mejor_actual){
			mejor_clase = i;
			mejor_actual = acumulado[i];
		}
	}
	/*cout << mejor_clase << endl;
	if (mejor_clase != 0){
		char ch;
		cin >> ch;
	}*/


	return mejor_clase;
}

//-----------------------------------------------------------------------------------------------------
void Pattern::TestearRecursivoVariosDicionarios(const example_set &Es, const VectorVar &V, const Pattern &P3, const VectorVar &V3, const Pattern &P2, const VectorVar &V2, TestResult & result){

  /*int n_clases = V.SizeDomain(V.Consecuente());
  result.cubiertos.clear();
  result.no_cubiertos.clear();
  result.acc.clear();
  result.error_intrinseco_porClase.clear();

  for (int i=0; i<n_clases; i++){
    result.cubiertos.push_back(0);
    result.no_cubiertos.push_back(0);
    result.acc.push_back(0);
    result.error_intrinseco_porClase.push_back(0);
  }

  double aciertos = 0, noCub = 0;
  cout << endl;
  int old_value=0, new_value;
  for (int i=0; i<Es.N_Examples(); i++){
    string regla="";
    double umbral = 0.3;
    int clase_predicha = -1;
    new_value = 100*i/Es.N_Examples();
    if (new_value != old_value){
      deleteLine2();
      cout << "Procesando Ejemplos " << 1.0*new_value << " %" <<endl;
      old_value = new_value;
    }
		string antecedente1,antecedente2,antecedente3;
    TestearRecursivoUnEjemplo(Es,V,i,regla,0,1,umbral,clase_predicha, false, antecedente1);
    //P2.TestearRecursivoUnEjemplo(Es,V2,i,regla,0,1,umbral,clase_predicha);
    if (clase_predicha == -1){
      umbral = 0.3;
      P3.TestearRecursivoUnEjemplo(Es,V3,i,regla,0,1,umbral,clase_predicha, false, antecedente2);
      if (clase_predicha == -1){
        umbral = 0.3;
        P2.TestearRecursivoUnEjemplo(Es,V2,i,regla,0,1,umbral,clase_predicha, false, antecedente3);
				if (clase_predicha == -1){
					clase_predicha = P2.TesteoDistanciaHamming(Es, V2, i);

				}
      }
    }
    //char ch;
    //cin >> ch;
    int clase_real = Es.Data(i,V.Consecuente());
    result.cubiertos[clase_real]++;
    if (clase_predicha==clase_real){
      result.acc[clase_real]++;
      aciertos++;
    }
    else if (clase_predicha == -1){
      result.no_cubiertos[clase_real]++;
      noCub++;
    }
  }
  cout << "Aciertos: " << 100*aciertos/Es.N_Examples() << endl;
  cout << "Aciertos no cubiertos: " << 100*aciertos/(Es.N_Examples()-noCub) << endl;
  cout << "Porcent no cubiertos: " << 100*noCub/Es.N_Examples() << endl;


  ProcesarResultados(result);*/
}


//-----------------------------------------------------------------------------------------------------

void Pattern::TestearPatronesBasicos(const example_set &Es, const VectorVar &V, TestResult & result){

  int n_clases = V.SizeDomain(V.Consecuente());
  result.cubiertos.clear();
  result.no_cubiertos.clear();
  result.acc.clear();
  result.error_intrinseco_porClase.clear();

  for (int i=0; i<n_clases; i++){
    result.cubiertos.push_back(0);
    result.no_cubiertos.push_back(0);
    result.acc.push_back(0);
    result.error_intrinseco_porClase.push_back(0);
  }

  for (int i=0; i<Es.N_Examples(); i++){
    string aux="";
    for (int j=0; j<V.N_Antecedente(); j++){
      string auxaux="";
      for (int l=0; l<V.SizeDomain(j); l++){
        auxaux.push_back('0');
      }
      //auxaux.push_back(' ');
      auxaux[BetterEtiqueta(V,j,Es,i)]='1';
      aux += auxaux;
    }

    int clase = Es.Data(i,V.Consecuente());

    auto it = diccionario.find(aux);
    int mayor;
    if (it!=diccionario.end()){
      result.cubiertos[clase]++;
      mayor=0;
      for (int k=1; k<n_clases; k++){
        if (it->second.conseq[k]>it->second.conseq[mayor])
        mayor=k;
      }
      if (mayor==clase)
        result.acc[mayor]++;
    }
    else{
      result.no_cubiertos[clase]++;
    }
  }

  ProcesarResultados(result);
}

//-----------------------------------------------------------------------------------------------------

int  Pattern::TestearUnEjemploPatronesBasicos(const example_set &Es, int eje, const VectorVar &V, double &peso){

  int n_clases = V.SizeDomain(V.Consecuente());
  int result;
  peso = 0;

  string aux; aux.reserve(10); aux ="";
  for (int j=0; j<V.N_Antecedente(); j++){
    string auxaux;
		auxaux.reserve(V.SizeDomain(j));
    for (int l=0; l<V.SizeDomain(j); l++){
      auxaux.push_back('0');
    }
    //auxaux.push_back(' ');
    auxaux[BetterEtiqueta(V,j,Es, eje)]='1';
    aux += auxaux;
  }


  //cout << aux << endl;
  int clase = Es.Data(eje,V.Consecuente());

  auto it = diccionario.find(aux);
  int mayor;
  if (it!=diccionario.end()){
    double ejemplos = it->second.conseq[0];
    mayor=0;
    for (int k=1; k<n_clases; k++){
      ejemplos += it->second.conseq[k];
      if (it->second.conseq[k]>it->second.conseq[mayor])
      mayor=k;
    }
    result = mayor;
    peso = 1.0 * (it->second.conseq[mayor]/ (ejemplos));
  }
  else{
    result = -1;
  }

	return result;

}


//-----------------------------------------------------------------------------------------------------
string PonerTiempo2 (double x){
	int t = x;
	int horas = t / 3600;
	t = t - horas * 3600;
	int minutos = t / 60;
	t = t - minutos * 60;
	int seg = t;
	string sseg, smin, shoras;
	if (seg<10){
	  sseg = "0"+to_string(seg);
  }
	else {
		sseg = to_string(seg);
	}

	if (minutos<10){
		smin = "0"+to_string(minutos);
	}
	else{
		smin = to_string(minutos);
	}


	string salida = to_string(horas) + ":" + smin + ":" + sseg;
	return salida;
}

//-----------------------------------------------------------------------------------------------------
void Pattern::TestearPatronesBasicosClassic(const example_set &Es, const VectorVar &V, TestResult & result){

  int n_clases = V.SizeDomain(V.Consecuente());
  result.cubiertos.clear();
  result.no_cubiertos.clear();
  result.acc.clear();
  result.error_intrinseco_porClase.clear();

  for (int i=0; i<n_clases; i++){
    result.cubiertos.push_back(0);
    result.no_cubiertos.push_back(0);
    result.acc.push_back(0);
    result.error_intrinseco_porClase.push_back(0);
  }

  // Calculo el peso y la clase asociada a cada regla
	vector<int> claseR;
	claseR.reserve(diccionario.size());

	vector<double> pesoR;
	pesoR.reserve(diccionario.size());

  // Version calculo de peso normal
	/*int auxC;
	double sum;
	for(auto it = diccionario.begin(); it != diccionario.end(); it++){
		auxC = 0;
		sum = 0;
		sum = it->second.conseq[0];
		for (int i=1; i<n_clases; i++){
			sum += it->second.conseq[i];
			if (it->second.conseq[i]>it->second.conseq[auxC]){
				auxC=i;
			}
		}
		claseR.push_back(auxC);
		pesoR.push_back(it->second.conseq[auxC]/sum);
  }*/

	// Version calculo de peso MapReduce
	int auxC;
	double sum;
	for(auto it = diccionario.begin(); it != diccionario.end(); it++){
		auxC = 0;
		sum = 0;
		sum = it->second.conseq[0];
		for (int i=1; i<n_clases; i++){
			sum += it->second.conseq[i];
		}
		for (int i=0; i<n_clases; i++){
			if (2*(it->second.conseq[i]) - sum > 0){
				auxC=i;
			}
		}
		claseR.push_back(auxC);
		pesoR.push_back((it->second.conseq[auxC]-(sum-it->second.conseq[auxC]))/sum);
  }


  cout << "Calculadas clases y pesos\n\n";
  int cl;
	double percentTestSet = 0.001;  //Porcentaje sobre el total del test sobre el que se está experimentando
	int nex = Es.N_Examples() * percentTestSet;
	clock_t inicio = clock();
  for (int i=0; i<nex; i++){
		clock_t step = clock();
		long int porcion = (1.0*i/nex)*10000;
		double portion = porcion/100.0;
		clock_t projeccion=0;
		double t_uno = 0 ;
		if (i > 0){
			t_uno = 1.0* (step-inicio)/i;
			projeccion = t_uno * (nex-i);
		}

		deleteLine2();
		cout << "\tEx " << portion
				 << "\%  ETA: " << PonerTiempo2(projeccion/CLOCKS_PER_SEC)
				 << "   Uno: " << t_uno / CLOCKS_PER_SEC
				 << "   TiempoTotal: " <<  PonerTiempo2(t_uno * nex/CLOCKS_PER_SEC) << endl;

    int clase = Es.Data(i,V.Consecuente());
		cl = -1;
		double mejorAdapt=0;

		int j=0;
		for(auto it = diccionario.begin(); it != diccionario.end(); it++, j++){
			//double adaptAnt = V.Adaptacion(Es.Data(i), it->first, mejorAdapt);
			double adaptAnt = V.AdaptacionTNormProduct(Es.Data(i), it->first, 0);
			adaptAnt = adaptAnt * pesoR[j];
			if (adaptAnt > mejorAdapt){
				mejorAdapt = adaptAnt;
				cl = claseR[j];
			}
		}

		if (cl == -1){
			result.no_cubiertos[clase]++;
		}
		else {
			result.cubiertos[clase]++;
			if (cl == clase){
				result.acc[clase]++;
		  }
		}
  }

  ProcesarResultados(result);
}


//-------------------------------------------------------------------------------------------------
void Pattern::TestearPatronesBasicosClassicDisparos(const example_set &Es, const VectorVar &V, TestResult & result, vector<pair<int,pair <string, double> > > &disparos){

  int n_clases = V.SizeDomain(V.Consecuente());
  result.cubiertos.clear();
  result.no_cubiertos.clear();
  result.acc.clear();
  result.error_intrinseco_porClase.clear();

  for (int i=0; i<n_clases; i++){
    result.cubiertos.push_back(0);
    result.no_cubiertos.push_back(0);
    result.acc.push_back(0);
    result.error_intrinseco_porClase.push_back(0);
  }

  // Calculo el peso y la clase asociada a cada regla
	vector<int> claseR;
	claseR.reserve(diccionario.size());

	//vector<pair<int,double> >numeroReglaR;
	//numeroReglaR.reserve(Es.N_Examples());


	vector<double> pesoR;
	pesoR.reserve(diccionario.size());

  // Version calculo de peso normal
	/*int auxC;
	double sum;
	for(auto it = diccionario.begin(); it != diccionario.end(); it++){
		auxC = 0;
		sum = 0;
		sum = it->second.conseq[0];
		for (int i=1; i<n_clases; i++){
			sum += it->second.conseq[i];
			if (it->second.conseq[i]>it->second.conseq[auxC]){
				auxC=i;
			}
		}
		claseR.push_back(auxC);
		pesoR.push_back(it->second.conseq[auxC]/sum);
  }*/

	// Version calculo de peso MapReduce
	int auxC;
	double sum;
	for(auto it = diccionario.begin(); it != diccionario.end(); it++){
		auxC = 0;
		sum = 0;
		sum = it->second.conseq[0];
		for (int i=1; i<n_clases; i++){
			sum += it->second.conseq[i];
			if (it->second.conseq[auxC] < it->second.conseq[i]){
				auxC = i;
			}
		}
		claseR.push_back(auxC);
		double peso;
		bool PCF = true;
		if (PCF){
			peso = (1.0 * it->second.conseq[auxC] - (sum-it->second.conseq[auxC])) / (sum);
		}
		else {
			peso = (1.0 * it->second.conseq[auxC])/ (sum);
		}
		pesoR.push_back(peso);
  }


  //cout << "Calculadas clases y pesos\n\n";
  int cl;
	int nex = Es.N_Examples();
	clock_t inicio = clock();
  for (int i=0; i<Es.N_Examples(); i++){
		clock_t step = clock();
		long int porcion = (1.0*i/nex)*10000;
		double portion = porcion/100.0;
		clock_t projeccion=0;
		double t_uno = 0 ;
		if (i > 0){
			t_uno = 1.0* (step-inicio)/i;
			projeccion = t_uno * (nex-i);
		}

		deleteLine2();
		cout << "\tEx " << portion
				 << "\%  ETA: " << PonerTiempo2(projeccion/CLOCKS_PER_SEC)
				 << "   Uno: " << t_uno / CLOCKS_PER_SEC
				 << "   TiempoTotal: " <<  PonerTiempo2(t_uno * nex/CLOCKS_PER_SEC) << endl;

    int clase = Es.Data(i,V.Consecuente());
		cl = -1;
		double mejorAdapt=0;

		int j=0;
		pair<int,pair<string, double> > elegida;
		elegida.first = -1;
		elegida.second.first = "-";
		elegida.second.second = 0;
		disparos.push_back(elegida);
		for(auto it = diccionario.begin(); it != diccionario.end(); it++, j++){
			//double adaptAnt = V.Adaptacion(Es.Data(i), it->first, mejorAdapt, pesoR[j]);
			double adaptAnt = V.AdaptacionTNormProduct(Es.Data(i), it->first, mejorAdapt, pesoR[j]); // Version CON poda
			//double adaptAnt = V.AdaptacionTNormProduct(Es.Data(i), it->first, 0); // Version SIN poda
			//adaptAnt = adaptAnt * pesoR[j]; // En la version sin poda no se combina con el peso, por eso esta linea hay que descomentarla en ese caso.
			if (adaptAnt > mejorAdapt){
				mejorAdapt = adaptAnt;
				cl = claseR[j];
				elegida.first = cl;
				elegida.second.first = (it->first);
				elegida.second.second = mejorAdapt;
				disparos[i]=elegida;
			}
		}

		if (cl == -1){
			result.no_cubiertos[clase]++;
		}
		else {
			result.cubiertos[clase]++;
			if (cl == clase){
				result.acc[clase]++;
		  }
		}
  }

  ProcesarResultados(result);
}


//-------------------------------------------------------------------------------------------------
int Pattern::TestearPatronesBasicosClassic_UnEjemplo(const example_set &Es, const VectorVar &V, const int eje, const vector<double> &pesoR, const vector<int> &claseR, double &mu, string &antecedente){
		int cl = -1;
		mu=0;

		int j=0;
		for(auto it = diccionario.begin(); it != diccionario.end(); it++, j++){
			//double adaptAnt = V.Adaptacion(Es.Data(i), it->first, mejorAdapt, pesoR[j]);
			double adaptAnt = V.AdaptacionTNormProduct(Es.Data(eje), it->first, mu, pesoR[j]); // Version CON poda
			//double adaptAnt = V.AdaptacionTNormProduct(Es.Data(eje), it->first, 0); // Version SIN poda
			//adaptAnt = adaptAnt * pesoR[j]; // En la version sin poda no se combina con el peso, por eso esta linea hay que descomentarla en ese caso.
			if (adaptAnt > mu){
				mu = adaptAnt;
				cl = claseR[j];
				antecedente = it->first;
			}
		}

  return cl;
}


void Pattern::CalcularPesoYClases(vector<double> &pesoR, vector<int> &claseR, bool PCF){
	claseR.reserve(diccionario.size());
	pesoR.reserve(diccionario.size());

  if (!PCF){
		// Version calculo de peso normal
		int auxC;
		double sum;
		for(auto it = diccionario.begin(); it != diccionario.end(); it++){
			auxC = 0;
			sum = 0;
			sum = it->second.conseq[0];
			for (int i=1; i<n_clases; i++){
				sum += it->second.conseq[i];
				if (it->second.conseq[i]>it->second.conseq[auxC]){
					auxC=i;
				}
			}
			claseR.push_back(auxC);
			pesoR.push_back(it->second.conseq[auxC]/sum);
		}
	}
	else {
		// Version calculo de peso PCF
		int auxC;
		double sum;
		for(auto it = diccionario.begin(); it != diccionario.end(); it++){
			auxC = 0;
			sum = 0;
			sum = it->second.conseq[0];
			for (int i=1; i<n_clases; i++){
				sum += it->second.conseq[i];
			}
			for (int i=0; i<n_clases; i++){
				if (2*(it->second.conseq[i]) - sum > 0){
					auxC=i;
				}
			}
			claseR.push_back(auxC);
			pesoR.push_back((it->second.conseq[auxC]-(sum-it->second.conseq[auxC]))/sum);
		}
	}

}



//-----------------------------------------------------------------------------------------------------
int Pattern::N_Patrones_paraClase (int clase) const {
  if (diccionario.size()==0){
    return 0;
  }

  auto it = diccionario.begin();
  int n_clases = it->second.conseq.size();

  if (clase<0 or clase>=n_clases){
    return 0;
  }
  else {
    double sumatotal=0;
    while (it!=diccionario.end()){
      double suma=0;
      int mayor=0;
      for (int i=0; i<it->second.conseq.size(); i++){
        suma += it->second.conseq[i];
        if (it->second.conseq[mayor]<it->second.conseq[i]){
          mayor = i;
        }
      }
      if (mayor == clase){
        sumatotal++;
      }
      it++;
    }
    return sumatotal;
  }
}

//-----------------------------------------------------------------------------------------------------
int Pattern::N_Ejemplos_paraClase (int clase) const {
  if (diccionario.size()==0){
    return 0;
  }

  auto it = diccionario.begin();
  int n_clases = it->second.eje_list.size();

  if (clase<0 or clase>=n_clases){
    return 0;
  }
  else {
    double sumatotal=0;
    while (it!=diccionario.end()){
      double suma=0;
      int mayor=0;
      sumatotal += it->second.eje_list[clase];
      it++;
    }
    return sumatotal;
  }
}
//-----------------------------------------------------------------------------------------------------
double CalPeso2(vector<double> &v, int &clase){
	double suma=0, peso;
	int mayor=0;
	for (int i=1; i<v.size(); i++){
		suma += v[i];
		if (v[i]>v[mayor]){
			mayor = i;
		}
	}
	if (suma>0){
		clase = mayor;
		peso = mayor / (suma);
	}
	else{
		clase = -1;
		peso = -1;
	}
	return peso;
}

//-----------------------------------------------------------------------------------------------------
int Pattern::InferirDiccionarioPatrones (vectordouble &ejemplo, const VectorVar &V){
	auto it = diccionario.begin();
	double adapt_aux, adaptacion = 0, peso;
	int clase = -1;

	while (it!=diccionario.end() and adaptacion < 1 ){
		adapt_aux = V.Adaptacion(ejemplo,it->first, adaptacion);
		if (adapt_aux>0){
			int cl;
			peso = CalPeso2(it->second.conseq, cl);
			//adapt_aux = peso * adapt_aux;
			if (adapt_aux > adaptacion){
				adaptacion = adapt_aux;
				clase = cl;
			}
		}
		it++;
	}
	return clase;
}

//-----------------------------------------------------------------------------------------------------
double Pattern::InferirDicionario (const example_set &E, const VectorVar &V){
  double aciertos=0, errores=0, no_cubiertos=0;
  cout << endl;
  for (int i=0; i<E.N_Examples(); i++){
    deleteLine2();
    cout << "Ejemplo " << i << "/" << E.N_Examples() << endl;

    vectordouble w = E.Data(i);
    int clase_orig = E.Data(i,V.Consecuente());
    int clase_predicha = InferirDiccionarioPatrones(w,V);
    if (clase_predicha == -1){
      no_cubiertos++;
    }
    else if (clase_orig == clase_predicha){
      aciertos++;
    }
    else {
      errores++;
    }
  }
  cout << "                 Aciertos: " << 100*aciertos/E.N_Examples()<< endl;
  cout << "     Porcent No Cubiertos: " << 100*no_cubiertos/E.N_Examples() << endl;
  cout << " Aciertos sobre Cubiertos: " << 100*aciertos/(E.N_Examples()-no_cubiertos);
  return aciertos/E.N_Examples();
}


//-----------------------------------------------------------------------------------------------------