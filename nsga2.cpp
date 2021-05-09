#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <values.h>
#include <float.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <algorithm>
#include <unistd.h>

/***
    AUTHOR: Marlon da Silva Rogério
    DATA: 08 MAIO 2021
    DISCIPLINA: Tópicos Especiais em Inteligência Computacional 
    DESCRICAO: Non-Dominated Sorting Genetic Algorithm II| Otimizacao multiobjetivo | Minimização
***/

long geracoes = 10; //100; //numero de geracoes
const int tam_pop = 100; //tamanho da populacao
const int dimensoes_obj=2; //numero de objetivos
const int dimensoes_var=dimensoes_obj+10-1; //numero de variaveis de decisao (valor padrao do DTLZ2)
double limiteInferior=0, limiteSuperior=1; //todas as solucoes devem estar entre 0 e 1 (padrao do DTLZ2)
int objAtual=0; //usada para ordenar as solucoes de acordo com cada objetivo de forma simples no calculo da crowding distance
struct Individuo {
    double x[dimensoes_var];
    double fx[dimensoes_obj];
    int rank=1;
    double cd=0;
};

Individuo populacao[tam_pop*2]; //tamanho do arquivo é igual à população a=m

//prototipacao das funcoes
void inicializacao();
void aptidao();
void atribuirRanks();
void recombinacao();
void cruzamento();
void mutacao();

//problemas de otimizacao
void calcularDTLZ2(double* x, double* fx); //calcula a funcao objetivo de acordo com o problema da literatura DTLZ2
void calcularDTLZ3(double* x, double* fx); //calcula a funcao objetivo de acordo com o problema da literatura DTLZ3

//funcoes opcionais para ajudar no processamento do algoritmo
bool compararDuasSolucoes(Individuo a, Individuo b);//Funcao necessaria para usar o std::sort - faz com que ordene pelo rank e desempate pelo cd
bool compararDuasSolucoesCD(Individuo a, Individuo b);//Funcao necessaria para usar o std::sort - faz com que ordene de acordo com o valor do objetivo que está na variável global objAtual
bool eDominada(Individuo a, int rankAtual); //verifica se um individuo e dominado ou nao em relacao aos individuos de um dado rank

int main(const int argc, const char* argv[]){
    srand (time(NULL)); //inicializa a semente aleatória com o tempo atual

    inicializacao(); //linha 3 do pseudocodigo. definimos valores iniciais para as soluções, pode ser aleatorio
    
    for(long g=0;g<geracoes;g++){//laco principal
        // std::cout << "[ GERACAO " << g << " ]";
        aptidao(); //linha 6 do pseudocodigo. calculo da aptidao, fitness ou valor objetivo
        //linha 7 do pseudocodigo. nossa populacao ja e unica, entao nao juntamos
        //linha 8 do pseudocodigo, nao faremos. retornaremos a populacao atual
        atribuirRanks(); //linha 9 do pseudocodigo. atribui os ranks para cada individuo da populacao (nao vamos separar em conjuntos diferentes)
        //linha 10 nao fazemos, no nosso caso so vamos sobrescrever as solucoes ruins
        recombinacao(); //linhas 11 a 17 junta as duas populacoes, faz isso reordenando o vetor populacao, já que só os primeiros tam_pop indivíduos são usados como pais no próximo passo
        cruzamento(); //aqui esta implementado o cruzamento com a etapa de selecao dos pais junto
        mutacao(); //mutacao, acontece em cada variavel de decisao de cada individuo com uma dada probabilidade
    }
        
        for(int i=0;i<tam_pop;i++){
            for(int j=0;j<dimensoes_obj;j++)
                printf("%.3f ", populacao[i].fx[j]);
                printf(" -- ");
                for(int j=0;j<dimensoes_var;j++){
                    printf("%.3f ", populacao[i].x[j]);
                }
            printf("\n");
        }
}

void inicializacao(){
    for(int i=0;i<tam_pop*2;i++){ //vamos inicializar todas as solucoes para facilitar nossa implementacao
        for(int j=0;j<dimensoes_var;j++){
            populacao[i].x[j]=rand()/(double)RAND_MAX;
        }
    }
}

void aptidao(){
    for(int i=0;i<tam_pop*2;i++){
        calcularDTLZ2(populacao[i].x, populacao[i].fx);
//         calcularDTLZ3(populacao[i].x, populacao[i].fx);
    }
}

//atribui um rank a cada individuo da populacao. o rank e a sua fronteira nao dominada. nao vamos separar as fronteiras em vetores diferentes, vamos so atribuir o rank de cada individuo
void atribuirRanks(){
    //determina quais sao as solucoes nao dominadas, elas recebem rank 1. ignora essas de rank 1 e encontra todas as solucoes nao dominadas, elas recebem rank 2. ignora essas de rank 2 e encontra todas as solucoes nao dominadas, elas recebem rank 3 .... ate que todas as solucoes sejam classificadas                     
    
    // Laço externo | Ind_A a ser comparado com todos os outros de seu rank    
    for (int ind_A=0; ind_A<tam_pop*2; ind_A++){
        bool alfa = false;
        // Laço interno | Ind_B a se comparado ao Ind_A
        for (int ind_B = ind_A+1; ind_B<tam_pop*2; ind_B++){
            // Verifica se estão do mesmo rank 
            if (populacao[ind_A].rank == populacao[ind_B].rank){
                // Percorre todos os objetivos
                for (int obj=0; obj<dimensoes_obj; obj++){
                    // Se A > B no objetivo "obj", então "alfa" = true 
                    if (populacao[ind_A].fx[obj] > populacao[ind_B].fx[obj]){
                        alfa = true;
                        populacao[ind_B].rank ++;
                    // Se B > A no objetivo "obj", então "alfa" = false | Nesse caso não há necessidade de verificar os próximos objetivos 
                    }else if (populacao[ind_B].fx[obj] > populacao[ind_A].fx[obj]){
                        alfa = false;
                        break;
                    }
                }
            }
        }
        //Se alfa permacecer "false" indica que o Fx de Ind_A foi menor em pelo menos uma ocasião se comparado aos individuos de seu rank
        if (alfa == false){
            populacao[ind_A].rank ++;
        }
        // std::cout << "[ CONTADOR " << ind_A << " ]" << "ind_A rank " << populacao[ind_A].rank  << "\n";
    }
}

void crowdingDistance(int rank){ //calculo da crowding distance para as solucoes em um dado rank
    int indiceInicioRank = 0;
    int indiceFimRank = 0;
    int contador = 0;
    //Percorre a populacao
    for (int i=0; i<tam_pop*2; i++){
        //Verifica se os individuos do laço estão no rank informado
        if (populacao[i].rank == rank){
            contador++;
            //Quando, e somenete quando o contador for == a 1 o indiceInicioRank recebe o indice do laço 
            if(contador==1){
                indiceInicioRank = i;
            }
            //Enquanto o laço for incrementado para individuos do rank a variável indiceFimRank o indice i para garantir que ao final seja armazenado o ultimo indice 
            indiceFimRank = i;
            //zera a crowding distance de todas as solucoes do rank
            populacao[i].cd = 0;
        }
    }
    //Percorre todos os objetivos
    for (int obj=0; obj<dimensoes_obj; obj++){
        //Possibilidade de valores entre o maior e o menor valor para o objetivo atual
        double range = (populacao[indiceFimRank].fx[obj] - populacao[indiceInicioRank].fx[obj]);
        //Ordena os individuos do rank pela objetivo do laço atual
        std::sort(populacao+indiceInicioRank, populacao+indiceFimRank+1, compararDuasSolucoesCD);
        //Primeiro e ultimo recebe um valor infinito;
        populacao[indiceInicioRank].cd = DBL_MAX;
        populacao[indiceFimRank].cd = DBL_MAX;
        //Percorre o intervalo de individuos do rank ente o 2° e penultimo individuo
        for (int ind=indiceInicioRank + 1; ind < indiceFimRank - 1; ind ++){
            // CD = CD + (IND_POSTERIOR - IND_ANTERIOR)/range
            populacao[ind].cd = populacao[ind].cd + (populacao[ind+1].fx[obj] - populacao[ind-1].fx[obj])/range;
        }
        objAtual = obj;
    }

    //ordena as solucoes por um dos objetivos. a primeira e a ultima recebem cd infinito
    //as outras recebem a diferenca entre o valor objetivo da solucao que esta antes dela e da solucao que esta depois dela
    //ordena todas as solucoes pelo outro objetivo
    //repete
    //ate acabarem os objetivos
    
    //essa linha de codigo pode ser util
}

void recombinacao(){
    std::sort(populacao, populacao+tam_pop*2, compararDuasSolucoes);//ordenar as solucoes de acordo com o rank
    for(int i=0;i<populacao[tam_pop*2-1].rank;i++)//estamos calculando cd de todos os ranks, nem precisa disso tudo, o custo computacional aumenta, mas e mais facil de implementar assim
        crowdingDistance(i);
    std::sort(populacao, populacao+tam_pop*2, compararDuasSolucoes);//ordenar as solucoes de acordo com o rank e o cd
}

int selecionarPai(){
    int f1=rand()%tam_pop;
    int f2=rand()%tam_pop;
        
    if(populacao[f1].rank < populacao[f2].rank)
        return f1;
    
    if(populacao[f1].rank > populacao[f2].rank)
        return f2;
    
    //chegamos aqui porque e igual
    if(populacao[f1].cd > populacao[f2].cd)
        return f1;
    else
        return f2;
}

void cruzamento(){
    Individuo pai1, pai2, filho1, filho2;
    for(int i=0;i<tam_pop;i+=2){
        
        memcpy(&pai1, &populacao[selecionarPai()], sizeof(Individuo));
        memcpy(&pai2, &populacao[selecionarPai()], sizeof(Individuo));
        
        //CRUZAMENTO SBX
        double y1,y2, betaq;
		double distributionIndex=30.0;//Parametro da distribuicao de valores
        for(int i=0;i<dimensoes_var;i++){
            double x1=pai1.x[i];
            double x2=pai2.x[i];
            if(rand()/(double)RAND_MAX <= 0.5){
                if(abs(x1 - x2) > 1.0e-14){//menor diferenca permitida entre valores
                    if(x1<x2){
                        y1=x1;
                        y2=x2;
                    }else{
                        y1=x2;
                        y2=x1;
                    }
                    double rnd=rand()/(double)RAND_MAX;
                    double beta = 1.0 + (2.0 * (y1 - limiteInferior) / (y2 - y1));
                    double alpha = 2.0 - pow(beta, -(distributionIndex + 1.0));
                    
                    if (rnd <= (1.0 / alpha)) {
                        betaq = pow(rnd * alpha, (1.0 / (distributionIndex + 1.0)));
                    } else {
                        betaq = pow(1.0 / (2.0 - rnd * alpha), 1.0 / (distributionIndex + 1.0));
                    }
                    double c1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));

                    beta = 1.0 + (2.0 * (limiteSuperior - y2) / (y2 - y1));
                    alpha = 2.0 - pow(beta, -(distributionIndex + 1.0));

                    if (rnd <= (1.0 / alpha)) {
                        betaq = pow((rnd * alpha), (1.0 / (distributionIndex + 1.0)));
                    } else {
                        betaq = pow(1.0 / (2.0 - rnd * alpha), 1.0 / (distributionIndex + 1.0));
                    }
                    double c2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));
                    
                    if (c1 < limiteInferior)
                        c1 = limiteInferior;

                    if (c2 < limiteInferior)
                        c2 = limiteInferior;

                    if (c1 > limiteSuperior)
                        c1 = limiteSuperior;

                    if (c2 > limiteSuperior)
                        c2 = limiteSuperior;
                    
                    if(rand()/(double)RAND_MAX <= 0.5){
                        filho1.x[i]=c2;
                        filho2.x[i]=c1;
                    }else{
                        filho1.x[i]=c1;
                        filho2.x[i]=c2;
                    }
                }else{
                    filho1.x[i]=x1;
                    filho2.x[i]=x2;
                }
            }else{
                filho1.x[i]=x2;
                filho2.x[i]=x1;
            }
        }
        
        memcpy(&populacao[i+tam_pop], &filho1, sizeof(Individuo));
        memcpy(&populacao[i+tam_pop+1], &filho2, sizeof(Individuo));
    }
}

void mutacao(){
    double probabilidade=1.0/dimensoes_var;
	double distributionIndex=30.0;
	double rnd, delta1, delta2, mut_pow, deltaq;
	double y, yl, yu, val, xy;
    
    //Mutacao polinomial
    for(int i=tam_pop;i<tam_pop*2;i++){
        for (int var=0; var < dimensoes_var; var++) {
            if ((rand()/(double)RAND_MAX) <= probabilidade){
                y      = populacao[i].x[var];
                yl     = limiteInferior;
                yu     = limiteSuperior;
                delta1 = (y-yl)/(yu-yl);
                delta2 = (yu-y)/(yu-yl);
                rnd=(rand()/(double)RAND_MAX);
                mut_pow = 1.0/(distributionIndex+1.0);
                if (rnd <= 0.5){
                    xy     = 1.0-delta1;
                    val    = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(distributionIndex+1.0)));
                    deltaq =  pow(val,mut_pow) - 1.0;
                }
                else{
                    xy = 1.0-delta2;
                    val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(distributionIndex+1.0)));
                    deltaq = 1.0 - (pow(val,mut_pow));
                }
                y = y + deltaq*(yu-yl);
                if (y<yl)
                    y = yl;
                if (y>yu)
                    y = yu;
                populacao[i].x[var]= y;
            }
        }
    }
}

//Funcao necessaria para usar o std::sort
bool compararDuasSolucoes(Individuo a, Individuo b){
    if(dimensoes_obj == 1)
        return a.fx[0] < b.fx[0];
    
    if(a.rank != b.rank)
        return a.rank < b.rank;
    
    return a.cd > b.cd; //CD quanto maior, melhor
}

//Funcao necessaria para usar o std::sort
bool compararDuasSolucoesCD(Individuo a, Individuo b){
    return a.fx[objAtual] < b.fx[objAtual];
}

bool eDominada(Individuo a, int rankAtual){
    //compara um individuo com todos os outros do mesmo rank e retorna se e dominado ou nao
    return false; //tem que ter um return pra compilar e nao dar problema
}

void calcularDTLZ2(double* x, double* fx){
	int k= dimensoes_var - dimensoes_obj + 1;
	double g=0.0;

	for(int i=dimensoes_var-k;i<dimensoes_var;i++)
		g+=(x[i]-0.5)*(x[i]-0.5);

	for(int i=0;i<dimensoes_obj;i++)
		fx[i] = (1.0+g);

	for(int i=0;i<dimensoes_obj;i++){
		for(int j=0;j<dimensoes_obj-(i+1);j++)
			fx[i] *= cos(x[j]*0.5*M_PI);
        if(i != 0){
            int aux = dimensoes_obj - (i+1);
            fx[i] *= sin(x[aux]*0.5*M_PI);
        }
	}
}

void calcularDTLZ3(double* x, double* fx){
	int k= dimensoes_var - dimensoes_obj + 1;
	double g=0.0;

	for(int i=dimensoes_var-k;i<dimensoes_var;i++)
		g+=(x[i]-0.5)*(x[i]-0.5)-cos(20.0*M_PI*(x[i]-0.5));

	g=100.0*(k+g);
	for(int i=0;i<dimensoes_obj;i++)
		fx[i] = 1.0+g;

	for(int i=0;i<dimensoes_obj;i++){
		for(int j=0;j<dimensoes_obj-(i+1);j++)
			fx[i] *= cos(x[j]*0.5*M_PI);
			if(i != 0){
				int aux = dimensoes_obj - (i+1);
				fx[i] *= sin(x[aux]*0.5*M_PI);
			}
	}
}