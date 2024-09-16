from numpy import conjugate, real_if_close


def ThreeDimIntegral(MyFunc, Ranges = (1,1,1), NumberOfPoints = 100):
    SumOfIntegral = 0
    XRange = Ranges[0]
    YRange = Ranges[1]
    ZRange = Ranges[2]
    x0, x1, y0, y1, z0, z1 = 0,0,0,0,0,0
    for dx in range(0, NumberOfPoints):
        x0 = ((2*XRange)/(NumberOfPoints+1)) * dx - XRange
        x1 = ((2*XRange)/(NumberOfPoints+1)) * (dx+1) - XRange
        
        for dy in range(0,NumberOfPoints):
            y0 = ((2*YRange)/(NumberOfPoints+1)) * dy - YRange
            y1 = ((2*YRange)/(NumberOfPoints+1)) * (dy+1) - YRange
            
            for dz in range(0,NumberOfPoints):
                z0 = ((2+ZRange)/(NumberOfPoints+1)) * dz - ZRange 
                z1 = ((2+ZRange)/(NumberOfPoints+1)) * (dz+1) - ZRange
                
                SumOfIntegral += ((z1 - z0) * (y1 - y0) * (x1 - x0)) * ((MyFunc(x0,y0,z0)+MyFunc(x1,y1,z1))/2)
                
    return SumOfIntegral

def Derivative(MyFunc, Coordinates, arg, d = 1e-7):
    x, y, z = Coordinates[0], Coordinates[1], Coordinates[2]
    if arg == 'x':
        DiffCoordinates = Coordinates.copy()
        DiffCoordinates[0]= Coordinates[0]+d
        return (MyFunc(DiffCoordinates) - MyFunc(Coordinates))/d
    elif arg == 'y':
        DiffCoordinates = Coordinates.copy()
        DiffCoordinates[1]= Coordinates[1]+d
        return (MyFunc(DiffCoordinates)-MyFunc(Coordinates))/d
    else:
        DiffCoordinates = Coordinates.copy()
        DiffCoordinates[2]= Coordinates[2]+d
        return (MyFunc(DiffCoordinates)-MyFunc(Coordinates))/d

def ModSquared(value):
    return real_if_close(value*conjugate(value))