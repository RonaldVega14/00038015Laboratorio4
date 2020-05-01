enum lines
{
    NOLINE,
    SINGLELINE,
    DOUBLELINE
};
enum modes
{
    NOMODE,
    INT_FLOAT,
    INT_INT_INT
};
enum parameters
{
    ELEMENT_LENGTH,
    PI,
    KAPPA,
    LAMBDA,
    IPSILON,
    PSI,
    ALPHA,
    DELTA,
    ETA
};
enum sizes
{
    NODES,
    ELEMENTS,
    DIRICHLET
};

class item
{

protected:
    int id;
    float x;
    int node1;
    int node2;
    float value;

public:
    int getId()
    {
        return id;
    }

    float getX()
    {
        return x;
    }

    int getNode1()
    {
        return node1;
    }

    int getNode2()
    {
        return node2;
    }

    float getValue()
    {
        return value;
    }

    void setId(int identifier)
    {
        id = identifier;
    }

    void setX(float x_coord)
    {
        x = x_coord;
    }

    void setNode1(int nodo)
    {
        node1 = nodo;
    }

    void setNode2(int nodo)
    {
        node2 = nodo;
    }

    void setValue(float node_val)
    {
        value = node_val;
    }

    virtual void setIntFloat(int n, float r) = 0;
    virtual void setIntIntInt(int n1, int n2, int n3) = 0;
};

class node : public item
{

public:
    void setIntFloat(int identifier, float x_coordinate)
    {
        id = identifier;
        x = x_coordinate;
    }

    void setIntIntInt(int n1, int n2, int n3)
    {
    }
};

class element : public item
{

public:
    void setIntFloat(int n1, float r)
    {
    }

    void setIntIntInt(int identifier, int firstnode, int secondnode)
    {
        id = identifier;
        node1 = firstnode;
        node2 = secondnode;
    }
};

class condition : public item
{

public:
    void setIntFloat(int node_to_apply, float prescribed_value)
    {
        node1 = node_to_apply;
        value = prescribed_value;
    }

    void setIntIntInt(int n1, int n2, int n3)
    {
    }
};

class mesh
{
    float parameters[9];
    int sizes[3];
    node *node_list;
    element *element_list;
    condition *dirichlet_list;

public:
    void setParameters(float l, float t, float k, float lamda,
                       float v, float psy, float alpha, float d,
                       float n)
    {
        parameters[ELEMENT_LENGTH] = l;
        parameters[PI] = t;
        parameters[KAPPA] = k;
        parameters[LAMBDA] = lamda;
        parameters[IPSILON] = v;
        parameters[PSI] = psy;
        parameters[ALPHA] = alpha;
        parameters[DELTA] = d;
        parameters[ETA] = n;
    }
    void setSizes(int nnodes, int neltos, int ndirich)
    {
        sizes[NODES] = nnodes;
        sizes[ELEMENTS] = neltos;
        sizes[DIRICHLET] = ndirich;
    }
    int getSize(int s)
    {
        return sizes[s];
    }
    float getParameter(int p)
    {
        return parameters[p];
    }
    void createData()
    {
        node_list = new node[sizes[NODES]];
        element_list = new element[sizes[ELEMENTS]];
        dirichlet_list = new condition[sizes[DIRICHLET]];
    }
    node *getNodes()
    {
        return node_list;
    }
    element *getElements()
    {
        return element_list;
    }
    condition *getDirichlet()
    {
        return dirichlet_list;
    }
    node getNode(int i)
    {
        return node_list[i];
    }
    element getElement(int i)
    {
        return element_list[i];
    }
    condition getCondition(int i, int type)
    {
        return dirichlet_list[i];
    }
};