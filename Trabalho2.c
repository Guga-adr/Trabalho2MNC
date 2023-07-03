// Autores
// Gustavo Amaral Duarte Rego
// Natan Mendes Alcantara
// Vinicius Person de Oliveira

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MAX 100

float determinante(int ordem, float matriz[MAX][MAX]) {
    float det = 0;
    float cofator[MAX][MAX];
    int sinal = 1;

    // Caso base para matriz de ordem 1
    if (ordem == 1) {
        return matriz[0][0];
    }

    // Loop para percorrer a primeira linha da matriz
    for (int i = 0; i < ordem; i++) {
        int submatrizLinha = 0;
        int submatrizColuna = 0;

        // Criacao da matriz cofatora
        for (int linha = 1; linha < ordem; linha++) {
            for (int coluna = 0; coluna < ordem; coluna++) {
                if (coluna != i) {
                    cofator[submatrizLinha][submatrizColuna] = matriz[linha][coluna];
                    submatrizColuna++;

                    // Incrementa a coluna da submatriz
                    if (submatrizColuna == ordem - 1) {
                        submatrizColuna = 0;
                        submatrizLinha++;
                    }
                }
            }
        }

        // Chamada recursiva para calcular o determinante da submatriz
        det += sinal * matriz[0][i] * determinante(ordem - 1, cofator);

        // Alterna o sinal para o próximo cofator
        sinal = -sinal;
    }

    return det;
}

int fatorial(int x) {
    if (x == 0)
        return 1;
    return (x * fatorial(x - 1));
}

float newtonGregory(int n, float tabela[][2], float x) {

    float h = fabs(tabela[0][0] - tabela[1][0]);

    for (int i = 1; i < n; i++) {

        if (fabs(tabela[i][0] - tabela[i - 1][0] - h) > 0.0001) {
            printf("Intervalos nao igualmente espacados!!!");
            return -1;
        }
    }

    float **delta = (float **)malloc(n * sizeof(float *));
    for (int i = 0; i < n; i++) {
        delta[i] = (float *)malloc(n * sizeof(float));
    }

    for (int i = 0; i < n; i++)
        delta[i][0] = tabela[i][1];

    for (int i = 1; i < n; i++) {

        for (int j = 0; j < n - i; j++) {

            delta[j][i] = (delta[j + 1][i - 1] - delta[j][i - 1]);
        }
    }

    float resultado = delta[0][0], termo = 1;

    for (int i = 1; i < n; i++) {

        for (int j = 0; j < i; j++)
            termo *= x - tabela[j][0];
        termo *= delta[0][i] / (fatorial(i) * pow(h, (double)i));
        resultado += termo;
        termo = 1;
    }

    for (int i = 0; i < n; i++) {
        free(delta[i]);
    }
    free(delta);

    return resultado;
}

float newton(int n, float tabela[][2], float x) {

    float **delta = (float **)malloc(n * sizeof(float *));
    for (int i = 0; i < n; i++) {
        delta[i] = (float *)malloc(n * sizeof(float));
    }

    for (int i = 0; i < n; i++)
        delta[i][0] = tabela[i][1];

    for (int i = 1; i < n; i++) {

        for (int j = 0; j < n - i; j++) {

            delta[j][i] = (delta[j + 1][i - 1] - delta[j][i - 1]) / (tabela[j + i][0] - tabela[j][0]);
        }
    }

    float resultado = delta[0][0], termo = 1;

    for (int i = 0; i < n - 1; i++) {

        for (int j = 0; j < i + 1; j++)
            termo *= x - tabela[j][0];
        termo *= delta[0][i + 1];

        resultado += termo;
        termo = 1;
    }

    for (int i = 0; i < n; i++) {
        free(delta[i]);
    }
    free(delta);

    return resultado;
}

float coefDet(int n, float tabela[][2], float yAjustados[]) {

    float yMedio = 0;

    for (int i = 0; i < n; i++) {
        yMedio += tabela[i][1];
    }

    yMedio /= n;

    float SQE = 0, SQT = 0;

    for (int i = 0; i < n; i++) {
        SQE += pow(yAjustados[i] - tabela[i][1], 2);
        SQT += pow(tabela[i][1] - yMedio, 2);
    }

    return 1 - SQE / SQT;
}

float ajusteReta(int n, float tabela[][2], float *a0, float *a1, float vetorY[], float *coefDeterminacao) {

    printf("Pontos tabelados: \n");
    for (int i = 0; i < n; i++) {
        printf("(%.4f, %.4f)\n", tabela[i][0], tabela[i][1]);
    }

    float somaX = 0, somaX2 = 0, somaY = 0, somaXY = 0;

    for (int i = 0; i < n; i++) {
        // Soma os valores de x, x^2, y e x*y
        somaX += tabela[i][0];
        somaX2 += pow(tabela[i][0], 2);
        somaY += tabela[i][1];
        somaXY += tabela[i][0] * tabela[i][1];
    }

    // Calcula o valor de a0 e a1

    *a1 = (n * somaXY - somaX * somaY) / (n * somaX2 - pow(somaX, 2));
    *a0 = (somaY - (*a1) * somaX) / n;

    // Calcula os valores de y ajustados

    for (int i = 0; i < n; i++) {
        vetorY[i] = *a0 + (*a1) * tabela[i][0];
    }

    // Calcula o coeficiente de determinacao

    *coefDeterminacao = coefDet(n, tabela, vetorY);

    return 0;
}

void SistemaTriangularInferior(int ordem, float coeficientes[MAX][MAX], float vetorIndInf[], float *vetorSolInf) {
    vetorSolInf[0] = vetorIndInf[0] / coeficientes[0][0];

    for (int i = 1; i < ordem; i++) {
        float soma = 0;
        for (int j = 0; j < i; j++) {
            soma += coeficientes[i][j] * vetorSolInf[j];
        }
        vetorSolInf[i] = (vetorIndInf[i] - soma) / coeficientes[i][i];
    }
}

void SistemaTriangularSuperior(int ordem, float coeficientes[MAX][MAX], float vetorIndSup[], float *vetorSolSup) {
    vetorSolSup[ordem - 1] = vetorIndSup[ordem - 1] / coeficientes[ordem - 1][ordem - 1];

    for (int i = ordem - 2; i >= 0; i--) {
        float soma = 0;
        for (int j = i + 1; j <= ordem - 1; j++) {
            soma += coeficientes[i][j] * vetorSolSup[j];
        }
        vetorSolSup[i] = (vetorIndSup[i] - soma) / coeficientes[i][i];
    }
}

int convergenciaAKMaior(int ordem, float matriz[MAX][MAX]) {
    for (int i = 1; i <= ordem; i++) {
        // teste

        if (determinante(i, matriz) <= 0) { // se devolver 0

            // teste:
            printf("\nNao converge\n");
            return 0; // nao converge
        }
    }
    // teste:
    printf("\nConverge!\n");
    return 1;
}

void auxCholesky(int ordem, float matriz[MAX][MAX], float matrizCholesky[MAX][MAX]) {

    for (int i = 0; i < ordem; i++) {
        for (int j = 0; j < ordem; j++) {
            // Calcular os elementos da Diagonal Principal
            if (i == j) {
                if (i == 0) {
                    matrizCholesky[i][j] = sqrt((double)matriz[i][j]);
                } else {
                    float soma = 0;

                    for (int k = 0; k < i; k++) {
                        soma += pow((double)matrizCholesky[i][k], 2);
                    }
                    matrizCholesky[i][j] = sqrt((double)matriz[i][j] - soma);
                }
            } else {
                // calcular os elementos fora da diagonal principal
                if (j < i) {
                    float soma = 0;
                    for (int k = 0; k < j; k++) {
                        // printf("\nsoma = matrizCholesky[%d][%d] * matrizCholesky[%d][%d]\n", i, k, j, k);
                        soma += matrizCholesky[i][k] * matrizCholesky[j][k];
                    }
                    // printf("\nmatrizCholesky[%d][%d] = (%f - %f) / %f\n", i, j, matriz[i][j], soma, matrizCholesky[j][j]);
                    matrizCholesky[i][j] = (matriz[i][j] - soma) / matrizCholesky[j][j];

                } else {
                    matrizCholesky[i][j] = 0;
                }
            }
        }
    }
}

void Cholesky(int ordem, float coeficientes[MAX][MAX], float vetorIndCholesky[], float *vetorSolCholesky) {
    float matrizL[MAX][MAX];

    // Verifica convergência do método
    if (convergenciaAKMaior(ordem, coeficientes) == 0) {
        printf("O sistema nao converge!\n");
        return;
    }
    auxCholesky(ordem, coeficientes, matrizL);
    printf("Matriz L: \n");
    for (int i = 0; i < ordem; i++) {
        for (int j = 0; j < ordem; j++) {
            printf("%.4f ", matrizL[i][j]);
        }
        printf("\n");
    }
    // Resolvendo o sistema Ly = b utilizando a funcao de matriz triangular inferior
    float vetorY[MAX];
    SistemaTriangularInferior(ordem, matrizL, vetorIndCholesky, vetorY);
    // Fazer a matriz transposta de L
    float matrizLT[MAX][MAX];
    for (int i = 0; i < ordem; i++) {
        for (int j = 0; j < ordem; j++) {
            matrizLT[i][j] = matrizL[j][i];
        }
    }
    // Resolvendo o sistema L^t x = y utilizando a funcao de matriz triangular superior
    SistemaTriangularSuperior(ordem, matrizLT, vetorY, vetorSolCholesky);
}

float ajustePolinomial(int n, int grauDesejado, float tabela[][2], float vetorA[], float vetorY[], float *coefDeterminacao) {

    // Inicia a Matriz de Minimos Quadrados

    float matrizMinimos[MAX][MAX];
    for (int i = 0; i < grauDesejado + 1; i++) {
        for (int j = 0; j < grauDesejado + 2; j++) {
            matrizMinimos[i][j] = 0;
        }
    }

    // calcular os somatórios necessários para a resolucao do sistema

    for (int i = 0; i < grauDesejado + 1; i++) {
        for (int j = 0; j < grauDesejado + 1; j++) {
            for (int k = 0; k < n; k++) {
                matrizMinimos[i][j] += pow(tabela[k][0], i + j);
            }
        }
    }

    // Monta um vetor Identidade como resolucao da matriz (somatórios X*Y) e atribiu-los ao vetorMinimos
    float vetorMinimos[MAX];
    for (int i = 0; i < n + 1; i++) {
        vetorMinimos[i] = 0;
        for (int j = 0; j < n + 1; j++) {
            vetorMinimos[i] += pow(tabela[j][0], i) * tabela[j][1];
        }
    }

    printf("\nMatriz de Minimos Quadrados:\n");

    for (int i = 0; i < grauDesejado + 1; i++) {
        for (int j = 0; j < grauDesejado + 1; j++) {
            printf("%.4f ", matrizMinimos[i][j]);
        }
        printf("\n");
    }

    printf("\nVetor de Minimos Quadrados:\n");

    for (int i = 0; i < grauDesejado + 1; i++) {
        printf("%.4f ", vetorMinimos[i]);
    }
    printf("\n");

    // Resolve o sistema de Minimos Quadrados por Cholesky e atribiu à vetorA

    Cholesky(grauDesejado + 1, matrizMinimos, vetorMinimos, vetorA);

    // Calcula o vetorY
    printf("Vetor A: \n");

    for (int i = 0; i < grauDesejado + 1; i++) {
        printf("%.4f ", vetorA[i]);
    }
    for (int i = 0; i < n; i++) {
        vetorY[i] = 0;
        for (int j = 0; j < grauDesejado + 1; j++) {
            vetorY[i] += vetorA[j] * pow(tabela[i][0], j);
        }
    }

    // Calcula o coeficiente de determinacao usando a funcao coefDeterminacao

    *coefDeterminacao = coefDet(n, tabela, vetorY);

    return 0;
}

void ajusteExponencial(int n, float tabela[][2], float *a, float *b, float vetorY[], float *coefDeterminacao) {
    // Inicia variáveis necessárias para o cálculo

    float somaX = 0, somaX2 = 0, somaY = 0, somaXY = 0;
    float a0 = 0, a1 = 0;

    // Calcular os somatórios necessários para a resolucao do sistema

    for (int i = 0; i < n; i++) {
        somaX += tabela[i][0];
        somaX2 += pow(tabela[i][0], 2);
        somaY += log(tabela[i][1]);
        somaXY += tabela[i][0] * log(tabela[i][1]);
    }

    // Calcula o valor de a0 e a1

    a1 = (n * somaXY - somaX * somaY) / (n * somaX2 - pow(somaX, 2));
    a0 = (somaY - a1 * somaX) / n;
    printf("a0 = %.4f\na1 = %.4f\n", a0, a1);
    *a = exp(a0);
    *b = exp(a1);
}

int main() {

    int exit = 0;
    int nPontosTabelados;
    float tabelaPontos[MAX][2], vetorY[MAX];
    do {

        printf("Digite o metodo de resolucao desejado: \n");
        printf("1 - Rotina Newton\n");
        printf("2 - Rotina Newton-Gregory\n");
        printf("3 - Rotina Coeficiente de Determinacao\n");
        printf("4 - Rotina Ajuste da Reta\n");
        printf("5 - Rotina Ajuste Polinomial\n");
        printf("6 - Rotina Ajuste Exponencial\n");
        printf("7 - Sair\n");

        int opcao;

        scanf("%d", &opcao);

        switch (opcao) {
            case 1:;
                nPontosTabelados = 0;
                printf("Digite o numero de pontos tabelados: ");
                scanf("%d", &nPontosTabelados);

                printf("Digite os pontos tabelados (X e depois Y)\n");
                for (int i = 0; i < nPontosTabelados; i++) {
                    for (int j = 0; j < 2; j++) {
                        scanf("%f", &tabelaPontos[i][j]);
                    }
                }
                float p;
                printf("Digite o valor de p: ");
                scanf("%f", &p);
                printf("Resultado: %.4f\n", newton(nPontosTabelados, tabelaPontos, p));

                break;
            case 2:;
                nPontosTabelados = 0;
                printf("Digite o numero de pontos tabelados: ");
                scanf("%d", &nPontosTabelados);

                printf("Digite os pontos tabelados (X e depois Y):\n");
                for (int i = 0; i < nPontosTabelados; i++) {
                    for (int j = 0; j < 2; j++) {
                        scanf("%f", &tabelaPontos[i][j]);
                    }
                }
                printf("Digite o valor de p: ");
                scanf("%f", &p);
                printf("Resultado: %.4f\n", newtonGregory(nPontosTabelados, tabelaPontos, p));

                break;
            case 3:;

                break;
            case 4:;
                nPontosTabelados = 0;
                float a0, a1;
                float coefDeterminacao = 0;
                printf("Digite o numero de pontos tabelados: ");
                scanf("%d", &nPontosTabelados);

                printf("Digite os pontos tabelados (X e depois Y):\n");
                for (int i = 0; i < nPontosTabelados; i++) {
                    for (int j = 0; j < 2; j++) {
                        scanf("%f", &tabelaPontos[i][j]);
                    }
                }
                ajusteReta(nPontosTabelados, tabelaPontos, &a0, &a1, vetorY, &coefDeterminacao);
                char sinal = '+';
                if (a1 < 0)
                    sinal = '-';
                printf("Equacao da reta: y = %.4f %c %.4f*X\n", a0, sinal, fabs(a1));
                printf("Coeficiente de determinacao: %.4f\n\n", coefDeterminacao);

                break;
            case 5:;
                nPontosTabelados = 0;
                int grauDesejado;
                printf("Digite o grau desejado: ");
                scanf("%d", &grauDesejado);

                printf("Digite o numero de pontos tabelados: ");
                scanf("%d", &nPontosTabelados);

                printf("Digite os pontos tabelados (X e depois Y):\n");
                for (int i = 0; i < nPontosTabelados; i++) {
                    for (int j = 0; j < 2; j++) {
                        scanf("%f", &tabelaPontos[i][j]);
                    }
                }

                float vetorA[MAX];
                ajustePolinomial(nPontosTabelados, grauDesejado, tabelaPontos, vetorA, vetorY, &coefDeterminacao);
                printf("\nVetor A: \n\n");
                for (int i = 0; i < grauDesejado + 1; i++) {
                    printf("%.4f ", vetorA[i]);
                }
                printf("\nEquacao do polinomio: ");
                for (int i = grauDesejado; i >= 0; i--) {
                    char sinal = '+';
                    if (vetorA[i] < 0)
                        sinal = '-';
                    if (i == 0)
                        printf("%.4f\n", vetorA[i]);
                    else
                        printf("%.4f*X^%d %c ", fabs(vetorA[i]), i, sinal);
                }

                printf("Vetor Y: \n\n");

                for (int i = 0; i < nPontosTabelados; i++) {
                    printf("%.4f ", vetorY[i]);
                }
                printf("\n");

                printf("\nCoeficiente de determinacao: %.4f\n\n", coefDeterminacao);

                break;
            case 6:;
                nPontosTabelados = 0;
                printf("Digite o numero de pontos tabelados: ");
                scanf("%d", &nPontosTabelados);

                printf("Digite os pontos tabelados (X e depois Y):\n");
                for (int i = 0; i < nPontosTabelados; i++) {
                    for (int j = 0; j < 2; j++) {
                        scanf("%f", &tabelaPontos[i][j]);
                    }
                }

                float a, b;
                ajusteExponencial(nPontosTabelados, tabelaPontos, &a, &b, vetorY, &coefDeterminacao);
                printf("Equacao exponencial: y = %.4f * %.4f^x\n", a, b);
                break;
            case 7:
                exit = 1;
                break;
            default:
                printf("Opcao invalida!\n");
                break;
        }
    } while (exit != 1);
}
