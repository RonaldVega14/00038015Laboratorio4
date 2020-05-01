void showVector(Vector b)
{
    cout << "[\t";
    for (int i = 0; i < b.size(); i++)
    {
        cout << b.at(i) << "\t";
    }
    cout << "]\n";
}

void showMatrix(Matrix K)
{
    for (int i = 0; i < K.at(0).size(); i++)
    {
        cout << "[\t";
        for (int j = 0; j < K.size(); j++)
        {
            cout << K.at(i).at(j) << "\t";
        }
        cout << "]\n";
    }
}

void showbs(vector<Vector> bs)
{
    for (int i = 0; i < bs.size(); i++)
    {
        cout << "b del elemento " << i + 1 << ":\n";
        showVector(bs.at(i));
        cout << "*************************************\n";
    }
}

void showKs(vector<Matrix> Ks)
{
    for (int i = 0; i < Ks.size(); i++)
    {
        cout << "K del elemento " << i + 1 << ":\n";
        showMatrix(Ks.at(i));
        cout << "*************************************\n";
    }
}