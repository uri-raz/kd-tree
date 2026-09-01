#define _CRT_SECURE_NO_WARNINGS

#include <iostream>

#include "kd tree.h"

int main()
{
	float x, y, z;
	std::vector<kd_tree<3>::kd_point> Points;

	auto fhandle = fopen("teapot_306.xyz", "r");
	while (fscanf(fhandle, "%f,%f,%f", &x, &y, &z) != EOF) {
		Points.push_back(kd_tree<3>::kd_point({ x, y, z }));
	}
	fclose(fhandle);

	kd_tree<3> Tree(Points);

	std::cout << "Tree has " << Tree.getNumPoints() << " points, and " << Tree.nodeCount() << " nodes.\n" << std::endl;

	kd_tree<3>::kd_point np = Tree.nearestNeighbor(kd_tree<3>::kd_point({ 3.0f, 0.0f, 2.5f })).value();

	std::cout << "Nearest point to (3, 0, 2.5) is (" << np[0] << ", " << np[1] << ", " << np[2] << ')' << std::endl;

	Points = Tree.KNearestNeighbors({ 3.0f, 0.0f, 2.5f }, 5).value();

	std::cout << "The 5nearest points to (3, 0, 2.5) are \n";
	for (const auto& p : Points) {
		std::cout << " (" << p[0] << ", " << p[1] << ", " << p[2] << ")\n";
	}
	std::cout << std::flush;

	return 0;
}
