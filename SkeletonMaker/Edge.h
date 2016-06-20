class Edge
{
public:
  Edge();
  ~Edge();

  float edgePoint[3];
  float a;
  float b;
  float w;
  int edgeTag;  // 0 for non-extremal, 1 for max, 2 for min
  int f;
  bool extremal, valid;

};
