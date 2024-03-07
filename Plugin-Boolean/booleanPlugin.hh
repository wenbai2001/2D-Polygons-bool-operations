#ifndef BOOLEANPLUGIN_HH_INCLUDED
#define BOOLEANPLUGIN_HH_INCLUDED
#include <OpenFlipper/BasePlugin/BaseInterface.hh>
#include <OpenFlipper/BasePlugin/ToolboxInterface.hh>
#include <OpenFlipper/BasePlugin/LoggingInterface.hh>
#include <OpenFlipper/BasePlugin/LoadSaveInterface.hh>
#include <OpenFlipper/common/Types.hh>

#include <ObjectTypes/PolyMesh/PolyMesh.hh>

#include <memory>
#include <tuple>
#include <QPushButton>
#include <QLabel>
#include <QGridLayout>
#include <QSpinBox>

class BooleanPlugin : public QObject, BaseInterface, ToolboxInterface, LoggingInterface, LoadSaveInterface
{
    Q_OBJECT
    Q_INTERFACES(BaseInterface)
    Q_INTERFACES(ToolboxInterface)
    Q_INTERFACES(LoggingInterface)
    Q_INTERFACES(LoadSaveInterface)
    Q_PLUGIN_METADATA(IID "org.OpenFlipper.Plugins.Boolean")

signals:
    void updateView();
    void updatedObject(int _identifier, const UpdateType &_type);

    void log(Logtype _type, QString _message);
    void log(QString _message);

    void addToolbox(QString _name, QWidget *_widget);

    void addEmptyObject(DataType _type, int &_id);

public:
    BooleanPlugin();

    QString name() { return QString("Boolean"); };

    QString description() { return QString("Does actually nothing but works!"); };

public:
    struct point {
        OpenMesh::Vec3d pos;
        bool isIntersection;
        bool in_out;
        bool is_visited;
    };
    using polygons = typename std::vector<point>;
    using polygons_ptr = typename std::shared_ptr<polygons>;

private:
    polygons_ptr Polygon1;
    polygons_ptr Polygon2;

    polygons_ptr load(const std::string& filePath);

    std::shared_ptr<PolyMesh> toMesh(polygons_ptr poly);

    void addMesh(std::shared_ptr<PolyMesh> polygon, const QString &name, int op);

private slots:
    void read_polygon1();
    void read_polygon2();
    void bool_intersection();
    void bool_union();
    void bool_difference_1();
    void bool_difference_2();
    void initializePlugin();

public slots:
    QString version() { return QString("Mesh-based"); };
};
#endif