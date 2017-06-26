#include <algorithm>
#include <cmath>
#include <iostream>
#include <queue>
#include <string>
#include <unordered_set>
#include <vector>

using namespace std;

typedef vector<double> vd;
typedef vector<int> vi;
typedef vector<vd> vvd;
typedef vector<char> vc;
typedef vector<bool> vb;


// Função hash para vectors, tirada do Stack Overflow
// Necessária para criar um unordered_set de vectors
class VectorHash {
public:
    size_t operator()(const vector<int>& v) const {
        hash<int> hasher;
        size_t seed = 0;
        for (int i : v) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        return seed;
    }
};

struct PPL{
    int n;                   // número de variáveis originais
    string tipo;             // "min" ou "max"
    string tipo_original;    // "min" ou "max"
    vd c;                    // vetor de custos
    vvd A;                   // matriz A
    vd b;                    // vetor b
    vi mod;                  // indica que tipo de modificação a variável mod i sofreu
                             //     0: variável inalterada
                             //     1: variável multiplicada por -1
                             //     2: variável trocada por xi - xi+1
    vc sinal_restricao;      // vetor com sinais de cada restrição
    vc sinal_variavel;       // vetor com sinais de cada variável
    vi folga_restricao;      // folga_restricao[i] é o índice da variável que
                             // é folga da restrição i, ou -1 se não foi adicionada folga.

    PPL() {}

    PPL(int n, int m) :
        n(n), c(n), A(m, vd(n)),  b(m), mod(n, 0),
        sinal_restricao(m), sinal_variavel(n), folga_restricao(m, -1)
        {}
};

struct Tabela_simplex {
    PPL& p;                  // referência para o PPL na forma padrão
    vi base;                 // vetor com os índices das variaveis básicas
    vi ibase;                // diz o índice da variável xi na base, ou -1 se xi não está na base
    vvd BiA;                 // matriz B^{-1}A
    vd Bib;                  // vetor B^{-1}b
    vd cB;                   // vetor c_B
    vd cz;                   // vetor de custos reduzidos
    double z;                // imagem atual, com o sinal trocado
    int primeira_artificial; // índice da primeira variável artificial
    int t;                   // tipo da solução (1..5). 0 se o algoritmo não terminou
    short fase;              // diz se estamos na fase 1 ou fase 2
    unordered_set<vi, VectorHash> otimos; // conjunto de vértices ótimos visitados
    queue<vi> fila_otimos;   // fila de vértices ótimos a serem visitados
    vvd solucao;             // vetor de soluções encontradas
    vvd direcao;             // direções extremas para laterais ótimas ilimitadas
    vi ivertice;             // diz o indíce do vértice associado a direção i

    Tabela_simplex(PPL& p) : 
        p(p), base(p.A.size()), ibase(p.A[0].size(), -1),
        Bib(p.b), cB(p.A.size()), t(0), fase(1), primeira_artificial(-1)
        {}
};

void ajustar_var_negativas(PPL &p);
void ajustar_var_irrestritas(PPL &p);
void ajustar_restricoes(PPL &p);
void ajustar_tipo(PPL &p);
void ajustar_tipo(PPL &p);
void forma_padrao(PPL &p);
void print_ppl(PPL &p, int k = 0);
bool atualizar(Tabela_simplex &tab);
int coluna_pivot(Tabela_simplex& tab);
int linha_pivot(int pivot, Tabela_simplex &tab);

bool deq(double a, double b) {
    return fabs(a-b) < 1e-7;
}

bool dgt(double a, double b) {
    return a > b and fabs(a-b) > 1e-7;
}

bool dge(double a, double b) {
    return dgt(a,b) or deq(a,b);
}

bool dlt(double a, double b) {
    return !dge(a,b);
}

void ajustar_var_negativas(PPL &p) {
    for (int i = 0; i < p.A[0].size(); ++i)
        if (p.sinal_variavel[i] == '<') {
            p.sinal_variavel[i] = '>';

            p.c[i] = -p.c[i];
            p.mod[i] = 1;

            for (int j = 0; j < p.A.size(); ++j)
                p.A[j][i] = -p.A[j][i];
        }
}

void ajustar_var_irrestritas(PPL &p) {
    for (int i = 0; i < p.A[0].size(); ++i)
        if (p.sinal_variavel[i] == 'L') {
            p.c.insert(p.c.begin() + i + 1, -p.c[i]);

            p.mod[i] = 2;
            for (int j = 0; j < p.A.size(); ++j)
                p.A[j].insert(p.A[j].begin() + i + 1, -p.A[j][i]);

            p.sinal_variavel[i] = '>';
            p.sinal_variavel.insert(p.sinal_variavel.begin() + i + 1, '>');
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

    string prefix("s.a ");
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

// ci - zi = ci - c_B^T B^{-1} a_i
double custo_reduzido(int j, Tabela_simplex &tab) {
    PPL& p = tab.p;
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

void mostrar_tabela(Tabela_simplex& tab) {
    cout.precision(4);
    cout << fixed;

    for (size_t i = 0; i < tab.cz.size(); ++i)
        cout << ((fabs(tab.cz[i]) < 5e-7)? 0.0: tab.cz[i]) << " ";
    cout << ((fabs(tab.z) < 5e-7)? 0.0: tab.z) << endl;

    for (size_t i = 0; i < tab.BiA.size(); ++i) {
        for (int j = 0; j < tab.BiA[0].size(); ++j)
            cout << ((fabs(tab.BiA[i][j]) < 5e-7)? 0.0: tab.BiA[i][j]) << " ";
        cout << ((fabs(tab.Bib[i]) < 5e-7)? 0.0: tab.Bib[i]) << '\n';
    }
}

void mostrar_solucao(Tabela_simplex& tab) {
    if (tab.t == 0) return;

    PPL &p = tab.p;

    cout << tab.t << endl;
    switch(tab.t) {
        case 1:
        case 2:
        case 3:
                for (size_t i = 0; i < tab.solucao.size(); ++i) {
                    cout << "V" << i+1;
                    for (size_t j = 0; j < tab.solucao[0].size(); ++j)
                        cout << " " << ((fabs(tab.solucao[i][j]) < 5e-7)? 0.0: tab.solucao[i][j]);
                    cout << '\n';
                }
                for (size_t i = 0; i < tab.direcao.size(); ++i) {
                    cout << "D" << tab.ivertice[i]+1;
                    for (size_t j = 0; j < tab.direcao[0].size(); ++j)
                        cout << " " << ((fabs(tab.direcao[i][j]) < 5e-7)? 0.0: tab.direcao[i][j]);
                    cout << '\n';
                }
                cout << (p.tipo_original == "max" ? tab.z : -tab.z) << endl;
                break;
        case 4: 
                cout << (p.tipo_original == "max" ? "+ inf" : "- inf") << endl;
                break;
        case 5: 
                cout << "imp" << endl;
    }
}

bool base_sem_artificial(Tabela_simplex& tab) {
    if (tab.primeira_artificial == -1)
        return true;

    for (int i = tab.primeira_artificial; i < tab.BiA[0].size(); ++i)
        if (tab.ibase[i] != -1)
            return false;

    return true;
}

int coluna_pivot(Tabela_simplex& tab) {
    int col_pivot = -1;
    double cz_min = INFINITY;

    switch(tab.fase) {
        case 1:
            for (int i = 0; i < tab.cz.size(); ++i) {
                if (tab.ibase[i] == -1 and dlt(tab.cz[i], cz_min)) {
                    cz_min = tab.cz[i];
                    col_pivot = i;
                }
            }

            return col_pivot;
        case 2:
            // regra de Bland: escolhe a primeira não-básica
            // com custo reduzido menor que zero
            for (int i = 0; i < tab.cz.size(); ++i) {
                if (tab.ibase[i] == -1 and dlt(tab.cz[i], 0.0))
                    return i;
                if (tab.ibase[i] == -1 and dlt(tab.cz[i], cz_min)) {
                    cz_min = tab.cz[i];
                    col_pivot = i;
                }
            }

            // se não existe não-básica com custo reduzido < 0
            // tentar colocar na fila vértices ótimos (custo 0) vizinhos
            vi base_viz(tab.base);
            for (int i = 0; i < tab.cz.size(); ++i) {
                if (tab.ibase[i] == -1 and deq(tab.cz[i], 0.0)) {
                    int lin_pivot = linha_pivot(i, tab);
                    if (lin_pivot == -1) continue;

                    base_viz[lin_pivot] = i;
                    if (tab.otimos.find(base_viz) == tab.otimos.end()) {
                        tab.fila_otimos.push(base_viz);
                        tab.otimos.insert(base_viz);
                    }
                    col_pivot = col_pivot == -1 ? i : col_pivot;
                    base_viz[lin_pivot] = tab.base[lin_pivot];
                }
            }
        return col_pivot;
    }
}

int linha_pivot(int pivot, Tabela_simplex &tab) {
    double razao = INFINITY;
    int lin_pivot = -1;

    for (size_t i = 0; i < tab.BiA.size(); ++i) {
        if (dgt(tab.BiA[i][pivot], 0)) {
            if (dgt(razao, tab.Bib[i] / tab.BiA[i][pivot])) {
                razao = tab.Bib[i] / tab.BiA[i][pivot];
                lin_pivot = i;
            } else if (deq(razao, tab.Bib[i] / tab.BiA[i][pivot]) and tab.base[lin_pivot] < tab.base[i]) {
                lin_pivot = i;
            }
        }
    }

    return lin_pivot;
}

void pivotear(int lin_pivot, int col_pivot, Tabela_simplex& tab) {
    double elem_pivot = tab.BiA[lin_pivot][col_pivot];

    tab.cB[lin_pivot] = tab.p.c[col_pivot];

    // atualiza linha pivot
    for (size_t i = 0; i < tab.BiA[0].size(); ++i)
        tab.BiA[lin_pivot][i] /= elem_pivot;
    tab.Bib[lin_pivot] /= elem_pivot;

    // atualiza outras linhas
    for (size_t i = 0; i < tab.BiA.size(); ++i)
        if (i != lin_pivot) {
            double alfa = -tab.BiA[i][col_pivot];
            for (int j = 0; j < tab.BiA[0].size(); ++j)
                tab.BiA[i][j] += alfa * tab.BiA[lin_pivot][j];

            tab.Bib[i] += alfa * tab.Bib[lin_pivot];
        }

    // atualiza custos reduzidos e imagem (com sinal trocado)
    for (size_t i = 0; i < tab.BiA[0].size(); ++i)
        tab.cz[i] = custo_reduzido(i, tab);
    tab.z = -imagem(tab);
}

void montar_solucao(Tabela_simplex& tab, int col_pivot = -2) {
    PPL& p = tab.p;
    int sidx = tab.solucao.size();

    cout.precision(4);
    cout << fixed;

    tab.otimos.insert(tab.base);

    tab.solucao.push_back(vd());
    for (int i = 0, o = 0; i < p.n; ++i, ++o) {
        int base_idx = tab.ibase[o];
        if (p.mod[i] == 1) {
            tab.solucao[sidx].push_back(base_idx != -1 ? -tab.Bib[base_idx] : 0.0);
        } else if (p.mod[i] == 2) {
            int base_idx2 = tab.ibase[o+1];
            tab.solucao[sidx].push_back((base_idx != -1 ? tab.Bib[base_idx] : 0.0) - (base_idx2 != -1 ? -tab.Bib[base_idx2] : 0.0));
            o++;
        } else {
            tab.solucao[sidx].push_back(base_idx != -1 ? tab.Bib[base_idx] : 0.0);
        }
    }

    if (col_pivot == -2) return;
    if (col_pivot == -1 and tab.t != 3) {
        tab.t = 2;
        return;
    }
    
    tab.t = 3;
    int didx = tab.direcao.size();

    tab.direcao.push_back(vd());
    tab.ivertice.push_back(sidx);
    for (int i = 0, o = 0; i < p.n; ++i, ++o) {
        int base_idx = tab.ibase[o];
        if (o == col_pivot)
            tab.direcao[didx].push_back(1.0);
        else if (p.mod[i] == 1) {
            tab.direcao[didx].push_back(base_idx != -1 ? -tab.BiA[base_idx][col_pivot] : 0.0);
        } else if (p.mod[i] == 2) {
            int base_idx2 = tab.ibase[o+1];
            tab.direcao[didx].push_back((base_idx != -1 ? -tab.BiA[base_idx][col_pivot] : 0.0) - (base_idx2 != -1 ? tab.BiA[base_idx2][col_pivot] : 0.0));
            o++;
        } else {
            tab.direcao[didx].push_back(base_idx != -1 ? -tab.BiA[base_idx][col_pivot] : 0.0);
        }
    }
}

bool atualizar(Tabela_simplex& tab) {
    PPL &p = tab.p;
    int col_pivot = coluna_pivot(tab);
    int lin_pivot = col_pivot != -1 ? linha_pivot(col_pivot, tab) : -1;
    double cz_min = col_pivot != -1 ? tab.cz[col_pivot] : INFINITY;

    switch (tab.fase) {
        case 1:
            if (deq(tab.z, 0) and base_sem_artificial(tab)) {
                return false;
            }
            else if (tab.z != 0 and dge(cz_min, 0)) {
                tab.t = 5;
                return false;
            }
            break;
        case 2:
            if (dgt(cz_min, 0) and tab.solucao.size() == 0) {
                tab.t = 1;
                montar_solucao(tab);
                return false;
            } else if (deq(cz_min, 0)) {
                montar_solucao(tab, lin_pivot == -1 ? col_pivot : -1);
            } else if (dlt(cz_min, 0) and lin_pivot == -1) {
                tab.t = 4;
                return false;
            }
            break;
    }

    if ((tab.t == 2 or tab.t == 3) and not tab.fila_otimos.empty()) {
        for (int i = 0; i < tab.base.size(); ++i)
            tab.ibase[tab.base[i]] = -1;

        tab.base = tab.fila_otimos.front();
        tab.fila_otimos.pop();

        for (int i = 0; i < tab.base.size(); ++i)
            tab.ibase[tab.base[i]] = i;

        for (int i = 0; i < tab.base.size(); ++i)
            pivotear(i, tab.base[i], tab);

        return true;
    } else if (not tab.t) {
        // muda a base
        tab.ibase[col_pivot] = lin_pivot;
        tab.ibase[tab.base[lin_pivot]] = -1;
        tab.base[lin_pivot] = col_pivot;

        pivotear(lin_pivot, col_pivot, tab);

        return true;
    } else {
        return false;
    }
}

void adicionar_var_artificiais(Tabela_simplex &tab) {
    PPL& p = tab.p;
    for (int i = 0; i < p.A.size(); ++i)
        // se a restrição i não tem folga, ou se a folga é negativa
        if (p.folga_restricao[i] == -1 or p.A[i][p.folga_restricao[i]] == -1) {
            if (tab.primeira_artificial == -1)
                tab.primeira_artificial = p.A[0].size();

            tab.base[i] = p.A[0].size();
            tab.ibase.push_back(i);

            for (int j = 0; j < p.A.size(); ++j)
                p.A[j].push_back(i == j ? 1 : 0);

            p.c.push_back(1);
            p.sinal_variavel.push_back('>');
        } else {
            tab.base[i] = p.folga_restricao[i];
            tab.ibase[p.folga_restricao[i]] = i;
        }
}

void remover_var_artificiais(Tabela_simplex& tab) {
    PPL& p = tab.p;
    if (tab.cz.size() > p.c.size()) {
        tab.cz.resize(p.c.size());
        for (int i = 0; i < tab.BiA.size(); ++i) {
            tab.BiA[i].resize(p.c.size());
        }
    }
}

void inicializar_tabela(Tabela_simplex &tab) {
    PPL& p = tab.p;
    if (tab.fase == 1)
        tab.BiA = p.A;

    for (int i = 0; i < p.A.size(); ++i)
        tab.cB[i] = p.c[tab.base[i]];

    if  (tab.cz.size() != p.c.size())
        tab.cz.resize(p.c.size());

    for (int i = 0; i < tab.cz.size(); i++)
        tab.cz[i] = custo_reduzido(i, tab);

    tab.z = -imagem(tab);
}

void fase1(Tabela_simplex &tab) {
    PPL &p = tab.p;
    adicionar_var_artificiais(tab);
    inicializar_tabela(tab);
    mostrar_tabela(tab);

    while (atualizar(tab));
}

void fase2(Tabela_simplex &tab) {
    if (tab.t)
        return;

    tab.fase = 2;
    remover_var_artificiais(tab);
    inicializar_tabela(tab);
    mostrar_tabela(tab);

    while (atualizar(tab));
}

void simplex(PPL &p) {
    Tabela_simplex tab(p);
    vd c_bak(p.c);
    
    p.c = vd(p.c.size(), 0);
    fase1(tab);

    p.c = c_bak;
    fase2(tab);

    mostrar_solucao(tab);
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
        p.tipo_original = p.tipo;
        for (int i = 0; i < n; ++i)
            cin >> p.c[i];

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j)
                cin >> p.A[i][j];
            cin >> p.sinal_restricao[i] >> p.b[i];
        }

        for (int i = 0; i < n; ++i) {
            cin >> p.sinal_variavel[i];
            if (p.sinal_variavel[i] == '<')
                p.mod[i] = 1;
            else if (p.sinal_variavel[i] == 'L')
                p.mod[i] = 2;
        }

        forma_padrao(p);
        if (!k) cout << "- - -\n";
        simplex(p);
        cout << "- - -\n";
    }
}