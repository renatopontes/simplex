#include <bits/stdc++.h>

using namespace std;

int main() {
	vector<int>	v(3);
	v[0] = 1;
	v[1] = 5;
	v[2] = 10;

	v.insert(v.begin() + 3, 99);

	// for (int i = 0; i < v.size(); ++i) {
	// 	cout << v[i] << ' ';
	// }
	cout << (v.end() - v.begin()) << endl;
}