#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

typedef vector<double> vd;
typedef vector<int> vi;
typedef vector<vd> vvd;
typedef vector<char> vc;
typedef vector<bool> vb;

struct PPL{
    string tipo;            // "min" ou "max"
    vd c;                   // vetor de custos
    vvd A;                  // matriz A
    vd b;                   // vetor b
    vc sinal_restricao;     // vetor com sinais de cada restrição
    vc sinal_variavel;      // vetor com sinais de cada variável
    vi folga_restricao;     // folga_restricao[i] é o índice da variável que é folga da restrição i, ou -1 se não foi adicionada folga.

    PPL() {}

    PPL(int n, int m) : c(n), A(m, vd(n)),  b(m), sinal_restricao(m), sinal_variavel(n), folga_restricao(m, -1) {}
};

struct Tabela_simplex {
    vi base;                // vetor com os índices das variaveis básicas
    vvd BiA;                // matriz B^{-1}A
    vd Bib;                 // vetor B^{-1}b
    vd cB;                  // vetor c_B
    vd cz;                  // vetor de custos reduzidos
    double z;               // imagem atual
    int solucao_encontrada;  // diz o tipo de solução da forma como a descrição do trabalho pede

    Tabela_simplex() {}

    Tabela_simplex(PPL& p) : base(p.A.size()), Bib(p.b), cB(p.A.size()), solucao_encontrada(0) {}
};

void ajustar_var_negativas(PPL &p);
void ajustar_var_irrestritas(PPL &p);
void ajustar_restricoes(PPL &p);
void ajustar_tipo(PPL &p);
void ajustar_tipo(PPL &p);
void forma_padrao(PPL &p);
void print_ppl(PPL &p, int k = 0);
void base_inicial(PPL &p, Tabela_simplex &tab);

void ajustar_var_negativas(PPL &p) {
    for (int i = 0; i < p.A[0].size(); ++i)
        if (p.sinal_variavel[i] == '<') {
            p.sinal_variavel[i] = '>';

            p.c[i] = -p.c[i];

            for (int j = 0; j < p.A.size(); ++j)
                p.A[j][i] = -p.A[j][i];
        }
}

void ajustar_var_irrestritas(PPL &p) {
    for (int i = 0; i < p.A[0].size(); ++i)
        if (p.sinal_variavel[i] == 'L') {
            p.c.insert(p.c.begin() + i + 1, -p.c[i]);

            for (int j = 0; j < p.A.size(); ++j)
                p.A[j].insert(p.A[j].begin() + i + 1, -p.A[j][i]);

            p.sinal_variavel.insert(p.sinal_variavel.begin() + i + 1, '>');
            p.sinal_variavel[i] = '>';
        }
}

void ajustar_restricoes(PPL &p) {
    for (int i = 0; i < p.A.size(); ++i)
        if (p.sinal_restricao[i] == '>' or p.sinal_restricao[i] == '<') {
            double folga = p.sinal_restricao[i] == '<' ? 1 : -1;

            p.folga_restricao[i] = p.A[0].size();

            for (int j = 0; j < p.A.size(); ++j)
                p.A[j].push_back(i == j ? folga : 0);

            p.c.push_back(0);
            p.sinal_variavel.push_back('>');

            p.sinal_restricao[i] = '=';
        }
}

void ajustar_tipo(PPL &p) {
    if (p.tipo == "max") {
        p.tipo = "min";

        for (int i = 0; i < p.A[0].size(); ++i)
            p.c[i] = -p.c[i];
    }
}

void forma_padrao(PPL &p) {
    ajustar_var_negativas(p);
    ajustar_var_irrestritas(p);
    ajustar_restricoes(p);
    ajustar_tipo(p);
}

void print_ppl(PPL &p, int k) {
    cout.precision(4);
    int n = p.A[0].size(), m = p.A.size();

    if (k) {
        cout << "P" << k << ":\n";
    }

    cout << p.tipo << ' ';
    for (int i = 0; i < n; ++i)
        cout << (i ? (p.c[i] >= 0.0 ? "+ " : "- ") : "") << (i ? fabs(p.c[i]) : p.c[i]) << "x" << i+1 << (i == n-1 ? "\n\n" : " ");

    string prefix("s.a. ");
    string fill(prefix.size(), ' ');

    cout << prefix;
    for (int i = 0; i < m; ++i) {
        cout << (i ? fill : "");
        for (int j = 0; j < n; ++j)
            cout << (j ? (p.A[i][j] >= 0.0 ? "+ " : "- ") : "") << (j ? fabs(p.A[i][j]) : p.A[i][j]) << "x" << j+1 << " ";
        cout << p.sinal_restricao[i] << (p.sinal_restricao[i] != '=' ? "=" : "") << " " << p.b[i] << (i == m-1 ? "\n\n" : "\n");
    }

    for (int i = 0; i < n; ++i) {
        cout << fill << "x" << i+1 << " ";
        if (p.sinal_variavel[i] == 'L')
            cout << "irrestrito\n";
        else
            cout << p.sinal_variavel[i] << "= 0\n";
    }
    cout << endl;
}

void adicionar_var_artificiais(PPL &p, Tabela_simplex &tab) {
    int m = p.A.size();
    for (int i = 0; i < m; ++i)
        if (p.folga_restricao[i] == -1 or p.A[i][p.folga_restricao[i]] == -1) {
            for (int j = 0; j < p.A.size(); ++j)
                p.A[j].push_back(i == j ? 1 : 0);

            p.c.push_back(1);
            p.sinal_variavel.push_back('>');
            tab.base[i] = p.A[0].size()-1;
        } else {
            tab.base[i] = p.folga_restricao[i];
        }
}

void remover_var_artificiais(PPL& p, Tabela_simplex& tab) {
    if (tab.cz.size() > p.c.size()) {
        tab.cz.resize(p.c.size());
        for (int i = 0; i < tab.BiA.size(); ++i) {
            tab.BiA[i].resize(p.c.size());
        }
    }
}


// ci - zi = ci - c_B^T B^{-1} a_i
double custo_reduzido(int j, Tabela_simplex &tab, PPL &p) {
    double r = 0;
    for (size_t i = 0; i < tab.cB.size(); ++i) {
        r += tab.cB[i] * tab.BiA[i][j];
    }

    return p.c[j] - r;
}

// z = c_B^T B^{-1} b
double imagem(Tabela_simplex &tab) {
    double r = 0;

    for (size_t i = 0; i < tab.cB.size(); ++i) {
        r += tab.cB[i] * tab.Bib[i];
    }

    return r;
}

void inicializar_tabela(PPL &p, Tabela_simplex &tab) {
    if (tab.BiA.size() == 0)
        tab.BiA = p.A;

    for (int i = 0; i < p.A.size(); ++i)
        tab.cB[i] = p.c[tab.base[i]];

    if  (tab.cz.size() != p.c.size())
        tab.cz.resize(p.c.size());

    for (int i = 0; i < tab.cz.size(); i++)
        tab.cz[i] = custo_reduzido(i, tab, p);

    tab.z = -imagem(tab);
}

void mostrar_tabela(Tabela_simplex &tab) {
    cout.precision(4);
    cout << fixed;

    for (size_t i = 0; i < tab.cz.size(); ++i)
        cout << tab.cz[i] << " ";
    cout << tab.z << endl;

    for (size_t i = 0; i < tab.BiA.size(); ++i) {
        for (int j = 0; j < tab.BiA[0].size(); ++j)
            cout << tab.BiA[i][j] << " ";
        cout << tab.Bib[i] << '\n';
    }
}

int teste_razao(int pivot, Tabela_simplex &tab) {
    double razao = INFINITY;
    int lin_pivot = -1;

    for (size_t i = 0; i < tab.BiA.size(); ++i) {
        if (tab.BiA[i][pivot] > 0 and tab.Bib[i] / tab.BiA[i][pivot] < razao) {
            razao = tab.Bib[i] / tab.BiA[i][pivot];
            lin_pivot = i;
        }
    }

    return lin_pivot;
}

void atualizar(Tabela_simplex &tab, PPL &p) {
    int col_pivot = min_element(tab.cz.begin(), tab.cz.end()) - tab.cz.begin();
    double cz_min = tab.cz[col_pivot];

    if (cz_min >= 0) {
        tab.solucao_encontrada = 1;
        return;
    }

    double lin_pivot = teste_razao(col_pivot, tab);
    if (lin_pivot == -1) return; // solução ilimitada

    double elem_pivot = tab.BiA[lin_pivot][col_pivot];

    tab.base[lin_pivot] = col_pivot;
    tab.cB[lin_pivot] = p.c[col_pivot];

    // atualiza linha pivot
    for (size_t i = 0; i < tab.BiA[0].size(); ++i)
        tab.BiA[lin_pivot][i] /= elem_pivot;
    tab.Bib[lin_pivot] /= elem_pivot;

    // atualiza outras linhas
    for (size_t i = 0; i < tab.BiA.size(); ++i)
        if (i != lin_pivot) {
            double alfa = -tab.BiA[i][col_pivot];
            for (int j = 0; j < tab.BiA[0].size(); ++j) {
                tab.BiA[i][j] += alfa * tab.BiA[lin_pivot][j];
            }
            tab.Bib[i] += alfa * tab.Bib[lin_pivot];
        }

    // atualiza custos reduzidos e imagem
    for (size_t i = 0; i < tab.BiA[0].size(); ++i)
        tab.cz[i] = custo_reduzido(i, tab, p);
    tab.z = imagem(tab);
}

void fase1(PPL &p, Tabela_simplex &tab) {
    adicionar_var_artificiais(p, tab);
    inicializar_tabela(p, tab);
    mostrar_tabela(tab);
    cout << endl;

    while (not tab.solucao_encontrada)
        atualizar(tab, p);
}

void fase2(PPL &p, Tabela_simplex &tab) {
    remover_var_artificiais(p, tab);
    inicializar_tabela(p, tab);
    mostrar_tabela(tab);

    while (not tab.solucao_encontrada)
        atualizar(tab, p);

    cout << "z: " << -tab.z << endl;
}

void simplex(PPL &p) {
    Tabela_simplex tab(p);
    vd c_bak(p.c);
    
    p.c = vd(p.c.size(), 0);
    fase1(p, tab);

    p.c = c_bak;
    fase2(p, tab);
}

int main() {
    int np;
    int indice = 1;

    cin >> np;

    for (int k = 0; k < np; ++k, ++indice) {
        int n, m;
        cin >> n >> m;

        PPL p(n, m);

        cin >> p.tipo;
        for (int i = 0; i < n; ++i)
            cin >> p.c[i];

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j)
                cin >> p.A[i][j];
            cin >> p.sinal_restricao[i] >> p.b[i];
        }

        for (int i = 0; i < n; ++i)
            cin >> p.sinal_variavel[i];

        forma_padrao(p);
        simplex(p);
    }
}