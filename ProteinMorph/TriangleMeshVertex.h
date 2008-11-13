#ifndef PROTEINMORPH_TRIANGLE_MESH_VERTEX_H
#define PROTEINMORPH_TRIANGLE_MESH_VERTEX_H

using namespace std;

namespace wustl_mm {
	namespace Protein_Morph {
		template <class TVertex> class TriangleMeshVertex {
		public:
			TriangleMeshVertex();
			TriangleMeshVertex(Vector3DFloat position);
		public:
			Vector3DFloat position;
			vector<unsigned long long> faceHashes;
			TVertex tag;
		};		

		template <class TVertex> TriangleMeshVertex<TVertex>::TriangleMeshVertex() {
			position = Vector3DFloat(0, 0, 0);
		}

		template <class TVertex> TriangleMeshVertex<TVertex>::TriangleMeshVertex(Vector3DFloat position) {
			this->position = position;
		}
	}
}

#endif