#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MAX 100

int main() {

    int exit = 0;
    do {
        double matriz[MAX][MAX];
        int ordem;
        printf("Digite a ordem da matriz: ");
        scanf("%d", &ordem);

        printf("Digite os elementos da matriz (Linha por linha):\n");
        for (int i = 0; i < ordem; i++) {
            for (int j = 0; j < ordem; j++) {
                scanf("%lf", &matriz[i][j]);
            }
        }

        printf("Digite o metodo de resolucao desejado: \n");
        printf("1 - Determinante\n");
        printf("2 - Triangular Inferior\n");
        printf("3 - Triangular Superior\n");
        printf("4 - Decomposicao LU\n");
        printf("5 - Rotina Cholesky\n");
        printf("6 - Rotina Gauss-Compacto\n");

        int opcao;

        scanf("%d", &opcao);

        switch (opcao) {
            case 1:
                printf("Determinante: %.4f\n", determinante(ordem, matriz));
                break;
            case 2:;
                double vetorIndInf[MAX];
                double vetorSolInf[MAX];
                printf("Digite os termos independentes: \n");
                for (int i = 0; i < ordem; i++) {
                    scanf("%lf", &vetorIndInf[i]);
                }
                SistemaTriangularInferior(ordem, matriz, vetorIndInf, vetorSolInf);
                printf("Vetor Solucao: \n");
                for (int i = 0; i < ordem; i++) {
                    printf("%.4f\n", vetorSolInf[i]);
                }
                break;
            case 3:;
                double vetorIndSup[MAX];
                double vetorSolSup[MAX];
                printf("Digite os termos independentes: \n");
                for (int i = 0; i < ordem; i++) {
                    scanf("%lf", &vetorIndSup[i]);
                }
                SistemaTriangularSuperior(ordem, matriz, vetorIndSup, vetorSolSup);
                printf("Vetor Solucao: \n");
                for (int i = 0; i < ordem; i++) {
                    printf("%.4f\n", vetorSolSup[i]);
                }
                break;
            case 4:;
                double vetorIndLU[MAX];
                double vetorSolLU[MAX];
                printf("Digite os termos independentes: \n");
                for (int i = 0; i < ordem; i++) {
                    scanf("%lf", &vetorIndLU[i]);
                }
                DecomposicaoLU(ordem, matriz, vetorIndLU, vetorSolLU);
                printf("Vetor Solucao: \n");
                for (int i = 0; i < ordem; i++) {
                    printf("%.4f\n", vetorSolLU[i]);
                }
                break;
            case 5:;
                double vetorIndCholesky[MAX];
                double vetorSolCholesky[MAX];
                printf("Digite os termos independentes: \n");
                for (int i = 0; i < ordem; i++) {
                    scanf("%lf", &vetorIndCholesky[i]);
                }
                Cholesky(ordem, matriz, vetorIndCholesky, vetorSolCholesky);
                printf("Vetor Solucao: \n");
                for (int i = 0; i < ordem; i++) {
                    printf("%.4f\n", vetorSolCholesky[i]);
                }
                break;
            case 6:;
                double vetorIndGaussCompacto[MAX];
                double vetorSolGaussCompacto[MAX];
                printf("Digite os termos independentes: \n");
                for (int i = 0; i < ordem; i++) {
                    scanf("%lf", &vetorIndGaussCompacto[i]);
                }
                GaussCompacto(ordem, matriz, vetorIndGaussCompacto, vetorSolGaussCompacto);
                printf("Vetor Solucao: \n");
                for (int i = 0; i < ordem; i++) {
                    printf("%.4f\n", vetorSolGaussCompacto[i]);
                }
                break;
            case 11:
                exit = 1;
                break;
            default:
                printf("Opcao invalida!\n");
                break;
        }
    } while (exit != 1);
}
