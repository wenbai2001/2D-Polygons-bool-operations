#include "booleanPlugin.hh"

#include <ObjectTypes/PolyMesh/PolyMesh.hh>

#include "OpenFlipper/BasePlugin/PluginFunctions.hh"

#include <iostream>
#include <cassert>
#include <string>
#include <fstream>
#include <vector>

// #define QUICKLOAD
int op_Polygon1 = -1;
int op_Polygon2 = -1;
int op_record = 0;
std::map<OpenMesh::Vec3d, bool> record1;
std::map<OpenMesh::Vec3d, bool> record2;

BooleanPlugin::BooleanPlugin() : Polygon1(nullptr), Polygon2(nullptr)
{

}

void BooleanPlugin::initializePlugin()
{
    QWidget *toolBox = new QWidget();

    QPushButton *loadPolygon1Btn = new QPushButton("&Polygon1", toolBox);
    QPushButton *loadPolygon2Btn = new QPushButton("&Polygon2", toolBox);
    QPushButton *EXECBooleanBtn = new QPushButton("&bool_intersection", toolBox);
    QPushButton *EXECBooleanBtn2 = new QPushButton("&bool_union", toolBox);
    QPushButton* EXECBooleanBtn3 = new QPushButton("&bool_difference_1", toolBox);
    QPushButton* EXECBooleanBtn4 = new QPushButton("&bool_difference_2", toolBox);

    QLabel *label = new QLabel("Iterations:");

    QGridLayout *layout = new QGridLayout(toolBox);
    layout->addWidget(loadPolygon1Btn, 0, 0);
    layout->addWidget(loadPolygon2Btn, 0, 1);
    layout->addWidget(EXECBooleanBtn, 1, 0, 1, 2);
    layout->addWidget(EXECBooleanBtn2, 2, 0, 1, 2);
    layout->addWidget(EXECBooleanBtn3, 3, 0, 1, 2);
    layout->addWidget(EXECBooleanBtn4, 4, 0, 1, 2);
    // layout->addItem(new QSpacerItem(10, 10, QSizePolicy::Expanding, QSizePolicy::Expanding), 2, 0, 1, 2);

    connect(loadPolygon1Btn, &QPushButton::clicked, this, &BooleanPlugin::read_polygon1);
    connect(loadPolygon2Btn, &QPushButton::clicked, this, &BooleanPlugin::read_polygon2);
    connect(EXECBooleanBtn, &QPushButton::clicked, this, &BooleanPlugin::bool_intersection);
    connect(EXECBooleanBtn2, &QPushButton::clicked, this, &BooleanPlugin::bool_union);
    connect(EXECBooleanBtn3, &QPushButton::clicked, this, &BooleanPlugin::bool_difference_1);
    connect(EXECBooleanBtn4, &QPushButton::clicked, this, &BooleanPlugin::bool_difference_2);
    emit addToolbox(tr("Boolean Operation"), toolBox);
}

void BooleanPlugin::read_polygon1()
{
#ifndef QUICKLOAD
    auto filePath = QFileDialog::getOpenFileName(nullptr,
                                                 tr("Choose one file"),
                                                 QString(),
                                                 u8"Plain Text (*.txt *.ppp)");
    Polygon1 = load(filePath.toStdString());
#else
    auto filePath = "res/Rect1.ppp";
    Polygon1 = load(filePath);
#endif
    auto mesh = toMesh(Polygon1);
    addMesh(mesh, tr("Polygon1"), 0);
}

void BooleanPlugin::read_polygon2()
{
#ifndef QUICKLOAD
    auto filePath = QFileDialog::getOpenFileName(nullptr,
                                                 tr("Choose one file"),
                                                 QString(),
                                                 u8"Plain Text (*.txt *.ppp)");
    Polygon2 = load(filePath.toStdString());
#else
    auto filePath = "res/Rect2.ppp";
    Polygon2 = load(filePath);
#endif
    auto mesh = toMesh(Polygon2);
    addMesh(mesh, tr("Polygon2"), 1);
}

// 检查两个浮点数是否足够接近
bool isCloseEnough(float a, float b) {
    return fabs(a - b) < 1e-6; // 1e-6是一个足够小的数，用来判断两个浮点数是否可以被认为是相等
}

//计算向量叉乘
float cross(const BooleanPlugin::point& O, const BooleanPlugin::point& A, const BooleanPlugin::point& B) {
    return (A.pos[0] - O.pos[0]) * (B.pos[1] - O.pos[1]) - (A.pos[1] - O.pos[1]) * (B.pos[0] - O.pos[0]);
}

//计算两点之间的距离
float distance(const BooleanPlugin::point& a, const BooleanPlugin::point& b) {
    return std::hypot(a.pos[0] - b.pos[0], a.pos[1] - b.pos[1]);
}

// 点是否在线段上
bool isPointOnSegment(BooleanPlugin::point P, BooleanPlugin::point Q, BooleanPlugin::point R) {
    if (Q.pos[0] <= std::max(P.pos[0], R.pos[0]) && Q.pos[0] >= std::min(P.pos[0], R.pos[0]) &&
        Q.pos[1] <= std::max(P.pos[1], R.pos[1]) && Q.pos[1] >= std::min(P.pos[1], R.pos[1])) {
        return true;
    }
    return false;
}

//插入交点到多边形中
void insertIntersection(std::vector<BooleanPlugin::point>& polygon, BooleanPlugin::point intersection) {
    // 为了插入交点，我们需要找到它应该插入的位置
    for (size_t i = 0; i < polygon.size(); ++i) {
        // 获取当前边的起点和终点
        BooleanPlugin::point start = polygon[i];
        BooleanPlugin::point end = polygon[(i + 1) % polygon.size()]; // 使用模运算来处理最后一个点到第一个点的情况

        // 检查交点是否在当前边上
        if (isCloseEnough(cross(start, intersection, end), 0.0) && isPointOnSegment(start, intersection, end)) {
            // 如果交点在当前边上，计算它与边起点的距离
            float distToStart = distance(start, intersection);
            float distToEnd = distance(end, intersection);

            // 在多边形顶点序列中找到正确的插入位置
            // 这个位置是基于交点与边起点之间的距离来确定的
            size_t insertPosition = i + 1;
            while (insertPosition < polygon.size() && distance(start, polygon[insertPosition]) < distToStart) {
                ++insertPosition;
            }

            // 将交点插入到顶点序列中的正确位置
            intersection.isIntersection = true; // 标记为交点
            polygon.insert(polygon.begin() + insertPosition, intersection);
            std::cout << "插入交点：" << intersection.pos << "\n";
            // 由于我们已经插入了交点，需要跳出循环
            break;
        }
    }
}

//求交点
bool intersectionPoint(OpenMesh::Vec3d A, OpenMesh::Vec3d B, OpenMesh::Vec3d C, OpenMesh::Vec3d D, OpenMesh::Vec3d& out) {
    double a1 = B[1] - A[1];
    double b1 = A[0] - B[0];
    double c1 = a1 * A[0] + b1 * A[1];

    double a2 = D[1] - C[1];
    double b2 = C[0] - D[0];
    double c2 = a2 * C[0] + b2 * C[1];

    double delta = a1 * b2 - a2 * b1;
    if (delta == 0) {
        // Lines are parallel
        return false;
    }

    out[0] = (b2 * c1 - b1 * c2) / delta;
    out[1] = (a1 * c2 - a2 * c1) / delta;
    out[2] = 0;
    if (out[0] < std::min(A[0], B[0]) || out[0] > std::max(A[0], B[0]) ||
        out[0] < std::min(C[0], D[0]) || out[0] > std::max(C[0], D[0]) ||
        out[1] < std::min(A[1], B[1]) || out[1] > std::max(A[1], B[1]) ||
        out[1] < std::min(C[1], D[1]) || out[1] > std::max(C[1], D[1])) {
        return false;
    }
    return true;
}

//点是否在多边形内，是1否0
bool isInside(const BooleanPlugin::point& test, const std::vector<BooleanPlugin::point>& polygon) {
    if (polygon.size() < 3) return false; // 不是多边形

    // 构造一个足够远的点，作为射线的终点
    std::vector<OpenMesh::Vec3d> extreme;
    OpenMesh::Vec3d extreme1 = { INT_MAX, test.pos[1], 0 };
    OpenMesh::Vec3d extreme2 = { INT_MAX, test.pos[1] + 10000, 0 };
    OpenMesh::Vec3d extreme3 = { test.pos[0], INT_MAX, 0 };
    OpenMesh::Vec3d extreme4 = { test.pos[0] - 10000, INT_MAX, 0 };
    extreme.push_back(extreme1);
    extreme.push_back(extreme2);
    extreme.push_back(extreme3);
    extreme.push_back(extreme4);
    OpenMesh::Vec3d inter;

    int in = 0;
    int out = 0;
    // 计数器，记录射线与多边形边的交点数
    for (int k = 0; k < 4; k++)
    {
        int count = 0, i = 0;
        do {
            int next = (i + 1) % polygon.size();

            // 检查多边形的边（polygon[i]，polygon[next]）是否与射线（test，extreme）相交
            if (intersectionPoint(polygon[i].pos, polygon[next].pos, test.pos, extreme[k], inter)) {
                count++;
            }
            std::cout << "check: "<<polygon[i].pos << " " << polygon[next].pos << " " << test.pos << " " << extreme[k] << " " << count << "\n";
            i = next;
        } while (i != 0);
        if (count % 2) in++;
        else out++;
    }
    if (in > out) return 1;
    return 0;
}

//标记进出点
void markEntryExitPoints(std::vector<BooleanPlugin::point>& poly1, std::vector<BooleanPlugin::point>& poly2, std::vector<BooleanPlugin::point> poly3, std::vector<BooleanPlugin::point> poly4)
{
    for (size_t i = 0; i < poly1.size(); ++i) {
        BooleanPlugin::point& p = poly1[i];
        if (p.isIntersection) {
            BooleanPlugin::point& prev = poly1[(i - 1 + poly1.size()) % poly1.size()];
            std::cout << "Poly1: Prev " << prev.pos[0] << prev.pos[1];
            if (prev.isIntersection) p.in_out = !prev.in_out;
            else
            {
                p.in_out = isInside(prev, poly4);
                std::cout << "not intersection " << p.in_out << "\n";
            }
        }
    }
    for (size_t i = 0; i < poly2.size(); ++i) {
        BooleanPlugin::point& p = poly2[i];
        if (p.isIntersection) {
            BooleanPlugin::point& prev = poly2[(i - 1 + poly2.size()) % poly2.size()];
            std::cout << "Poly2: Prev " << prev.pos[0] << prev.pos[1];
            if (prev.isIntersection) p.in_out = !prev.in_out;
            else
            {
                p.in_out = isInside(prev, poly3);
                std::cout << "not intersection " << p.in_out<<"\n";
            }
        }
    }
}

//多边形切换
BooleanPlugin::point* switchPolygon(BooleanPlugin::point* currentPoint, BooleanPlugin::polygons_ptr& poly1, BooleanPlugin::polygons_ptr& poly2, bool cur) {
    // 在另一个多边形中找到匹配的交点
    if (!cur)  //如果cur为0，当前为poly1，应该转到poly2上
    {
        for (auto& pt : *poly2) {
            if (pt.pos == currentPoint->pos) {
                return &pt;
            }
        }
    }
    else
    {
        for (auto& pt : *poly1) {
            if (pt.pos == currentPoint->pos) {
                return &pt;
            }
        }
    }
    return nullptr; // 如果没有找到匹配的交点，返回nullptr
}

//获取下一个
BooleanPlugin::point* getNextPoint(BooleanPlugin::point* currentPoint, BooleanPlugin::polygons_ptr& poly) {
    auto it = std::find_if(poly->begin(), poly->end(), [currentPoint](const BooleanPlugin::point& p) {
        return p.pos == currentPoint->pos;
        });
    if (it != poly->end()) {
        auto nextIt = std::next(it);
        if (nextIt == poly->end()) {
            nextIt = poly->begin();
        }
        return &(*nextIt);
    }
    return nullptr;
}

//获取上一个
BooleanPlugin::point* getPrevPoint(BooleanPlugin::point* currentPoint, BooleanPlugin::polygons_ptr& poly)
{
    for (size_t i = 0; i < poly->size(); ++i) {
        // 如果找到当前点
        if (&(*poly)[i] == currentPoint) {
            // 如果当前点是多边形的第一个点
            if (i == 0) {
                // 返回多边形的最后一个点的地址
                return &(*poly).back();
            }
            else {
                // 返回前一个点的地址
                return &(*poly)[i - 1];
            }
        }
    }
    // 如果没有找到当前点，返回空指针
    return nullptr;
}

//布尔交
void BooleanPlugin::bool_intersection()
{
    if (!Polygon1 || !Polygon2) {
        std::cout << "Load first" << std::endl;
        return;
    }

    std::vector<point> intersections;

    for (auto a_it = Polygon1->begin(); a_it != Polygon1->end(); ++a_it) {
        auto a_next = std::next(a_it);
        if (a_next == Polygon1->end()) a_next = Polygon1->begin();

        for (auto b_it = Polygon2->begin(); b_it != Polygon2->end(); ++b_it) {
            auto b_next = std::next(b_it);
            if (b_next == Polygon2->end()) b_next = Polygon2->begin();

            OpenMesh::Vec3d inter;
            if (intersectionPoint(a_it->pos, a_next->pos, b_it->pos, b_next->pos, inter)) {
                point p;
                p.pos = inter;
                p.isIntersection = 1;
                p.in_out = 0;
                p.is_visited = 0;
                intersections.push_back(p);
            }
        }
    }

    polygons polygon3, polygon4;
    polygon3 = *Polygon1;
    polygon4 = *Polygon2;


    for (int i = 0; i < intersections.size(); i++)
    {
        insertIntersection(*Polygon1, intersections[i]);
        insertIntersection(*Polygon2, intersections[i]);
    }
    markEntryExitPoints(*Polygon1, *Polygon2, polygon3, polygon4);
    
    if (Polygon1->begin()->pos[0] == -1.0 && Polygon1->begin()->pos[1] == 0.0)
    {
        return;
    }

    if (Polygon1->begin()->pos[0] == -1.0 && Polygon1->begin()->pos[1] == -1.0)
    {
        auto mesh = toMesh(Polygon2);
        addMesh(mesh, tr("Union2"), 2);
        return;
    }

    for (auto& point : *Polygon1) {
        polygons_ptr result = std::make_shared<polygons>();
        bool first = 0;
        bool cur_poly = 0;   //0表示当前在poly1上，应该转到2
        if (point.isIntersection && !point.in_out) { // 如果是进点
            auto currentPoint = &point;
            do {
                std::cout << currentPoint->pos[0] << " " << currentPoint->pos[1] << " " << currentPoint->in_out << "\n";
                result->emplace_back(*currentPoint);
                if(first)  //第一次进不需要转换多边形
                {
                    std::cout << "first in" << "\n";
                    if (currentPoint->isIntersection) {
                        currentPoint = switchPolygon(currentPoint, Polygon1, Polygon2, cur_poly);
                        cur_poly = !cur_poly;
                    }
                }
                if(cur_poly) currentPoint = getNextPoint(currentPoint, Polygon2);
                else currentPoint = getNextPoint(currentPoint, Polygon1);
                first = 1;
            } while (currentPoint->pos != point.pos); 
            auto mesh = toMesh(result);
            addMesh(mesh, tr("Intersection"), 1);
        }
    }

}

//布尔并
void BooleanPlugin::bool_union()
{
    if (Polygon1->begin()->pos[0] == -1.0 && Polygon1->begin()->pos[1] == 0.0)
    {
        polygons_ptr r = std::make_shared<polygons>();
        point po;
        auto mesh = toMesh(Polygon1);
        addMesh(mesh, tr("Union"), 2);
        mesh = toMesh(Polygon2);
        addMesh(mesh, tr("Union2"), 2);
        return;
    }
    if (Polygon1->begin()->pos[0] == -1.0 && Polygon1->begin()->pos[1] == -1.0)
    {
        auto mesh = toMesh(Polygon1);
        addMesh(mesh, tr("Union2"), 2);
        return;
    }
    polygons_ptr result = std::make_shared<polygons>();

    bool if_dec = 1;
    bool cur_poly = 0;   //0表示当前在poly1上，应该转到2
    for (auto& point : *Polygon1) {
        if (point.isIntersection && !point.in_out) { // 如果是进点
            if_dec = 0;
            auto currentPoint = &point;
            do {
                result->emplace_back(*currentPoint);

                if (currentPoint->isIntersection) {
                    currentPoint = switchPolygon(currentPoint, Polygon1, Polygon2, cur_poly);
                    cur_poly = !cur_poly;
                }
                if (cur_poly) currentPoint = getNextPoint(currentPoint, Polygon2);
                else currentPoint = getNextPoint(currentPoint, Polygon1);
            } while (currentPoint != &point);

            break;
        }
    }
    if (if_dec)
    {
        auto mesh = toMesh(Polygon1);
        addMesh(mesh, tr("Union"), 2);
        mesh = toMesh(Polygon2);
        addMesh(mesh, tr("Union2"), 2);
    }
    else
    {
        auto mesh = toMesh(result);
        addMesh(mesh, tr("Union"), 2);
    }
}

//布尔差1-2
void BooleanPlugin::bool_difference_1()
{
    if (Polygon1->begin()->pos[0] == -1.0 && Polygon1->begin()->pos[1] == 0.0)
    {
        auto mesh = toMesh(Polygon1);
        addMesh(mesh, tr("Union"), 2);
        return;
    }
    bool if_dec = 1;
    for (auto& point : *Polygon1) {
        polygons_ptr result = std::make_shared<polygons>();
        bool cur_poly = 0;   //0表示当前在poly1上，应该转到2
        if (point.isIntersection && !point.in_out) { // 如果是进点
            if_dec = 0;
            auto currentPoint = &point;
            do {
                std::cout << currentPoint->pos[0] << " " << currentPoint->pos[1] << " " << currentPoint->in_out << "\n";
                result->emplace_back(*currentPoint);
                if (currentPoint->isIntersection) {
                    currentPoint = switchPolygon(currentPoint, Polygon1, Polygon2, cur_poly);
                    cur_poly = !cur_poly;
                }
                if (cur_poly) currentPoint = getPrevPoint(currentPoint, Polygon2);
                else currentPoint = getNextPoint(currentPoint, Polygon1);
            } while (currentPoint->pos != point.pos);
            auto mesh = toMesh(result);
            addMesh(mesh, tr("Difference"), 1);
        }
    }
    if (if_dec)
    {
        auto mesh = toMesh(Polygon1);
        addMesh(mesh, tr("Difference"), 2);
    }
}

//布尔差2-1
void BooleanPlugin::bool_difference_2()
{
    if (Polygon1->begin()->pos[0] == -1.0 && Polygon1->begin()->pos[1] == 0.0)
    {
        auto mesh = toMesh(Polygon2);
        addMesh(mesh, tr("Union2"), 2);
        return;
    }
    if (Polygon1->begin()->pos[0] == -1.0 && Polygon1->begin()->pos[1] == -1.0)
    {
        return;
    }
    bool if_dec = 1;
    for (auto& point : *Polygon2) {
        polygons_ptr result = std::make_shared<polygons>();
        bool cur_poly = 0;   //0表示当前在poly2上，应该转到1
        if (point.isIntersection && !point.in_out) { // 如果是进点
            if_dec = 0;
            auto currentPoint = &point;
            do {
                std::cout << currentPoint->pos[0] << " " << currentPoint->pos[1] << " " << currentPoint->in_out << "\n";
                result->emplace_back(*currentPoint);
                if (currentPoint->isIntersection) {
                    cur_poly = !cur_poly;
                    currentPoint = switchPolygon(currentPoint, Polygon1, Polygon2, cur_poly);
                }
                if (cur_poly) currentPoint = getPrevPoint(currentPoint, Polygon1);  //顺时针遍历1
                else currentPoint = getNextPoint(currentPoint, Polygon2);
            } while (currentPoint->pos != point.pos);
            auto mesh = toMesh(result);
            addMesh(mesh, tr("Difference"), 1);
        }
    }
    if (if_dec)
    {
        auto mesh = toMesh(Polygon2);
        addMesh(mesh, tr("Difference"), 2);
    }
}

BooleanPlugin::polygons_ptr BooleanPlugin::load(const std::string &filePath)
{
    using namespace std;
    auto f = fstream(filePath, ios::in);
    auto ptr = make_shared<PolyMesh>();

    if (!f.is_open())
        return nullptr;

    string line, flag;
    istringstream linestream;

    std::vector<std::string> allLines;
    polygons_ptr r = std::make_shared<polygons>();

    // 读取并丢弃第一行
    std::getline(f, line);

    // 读取所有剩余行到临时容器
    while (std::getline(f, line)) {
        allLines.push_back(line);
    }

    // 关闭文件
    f.close();

    // 如果文件只有一行或没有行，则返回空的polygons_ptr
    if (allLines.size() <= 1) {
        return r;
    }

    // 跳过最后一行，将所有其他行转移到polygons数据结构中
    for (size_t i = 0; i < allLines.size() - 1; ++i) {
        std::istringstream linestream(allLines[i]);
        double x, y;
        const double z = 0; // 假设z坐标总是0
        linestream >> x >> y;
        point po;
        po.pos = OpenMesh::Vec3d(x, y, z);
        po.isIntersection = false;
        po.in_out = false;
        po.is_visited = false;
        r->emplace_back(po);
    }

    return r;
    /*
    polygons_ptr r = std::make_shared<polygons>();

    getline(f, line);
    while (getline(f, line))
    {
        linestream = istringstream(line);
        double x, y;
        const double z = 0;
        linestream >> x >> y;
        point po;
        po.pos = OpenMesh::Vec3d(x, y, z);
        po.isIntersection = 0;
        po.in_out = 0;
        po.is_visited = 0;
        r->emplace_back(po);
    }

    f.close();

    return r;*/

}

std::shared_ptr<PolyMesh> BooleanPlugin::toMesh(polygons_ptr poly)
{
    std::shared_ptr<PolyMesh> mesh = std::make_shared<PolyMesh>();

    std::vector<OpenMesh::SmartVertexHandle> v;
    for (auto it = poly->begin(); it != poly->end(); it++)
    {
        v.push_back(mesh->add_vertex(it->pos));
    }
    mesh->add_face(v);

    return mesh;
}

void BooleanPlugin::addMesh(std::shared_ptr<PolyMesh> polygon, const QString &name, int op)
{
    int newObjectId = -1;
    emit addEmptyObject(DATA_POLY_MESH, newObjectId);
    if (op) op_Polygon2 = newObjectId;
    else op_Polygon1 = newObjectId;
    PolyMeshObject *result = PluginFunctions::polyMeshObject(newObjectId);
    result->setName(name);
    result->mesh()->assign(*polygon, true);

    result->mesh()->request_vertex_status();
    result->mesh()->request_vertex_normals();
    result->mesh()->request_vertex_colors();
    result->mesh()->request_halfedge_status();
    result->mesh()->request_edge_status();
    result->mesh()->request_face_normals();
    result->mesh()->request_face_colors();
    result->mesh()->request_face_status();

    emit updatedObject(newObjectId, UPDATE_GEOMETRY);

}
