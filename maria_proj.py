from cmath import sqrt
import scipy.linalg as la



x,y,z = 30,10,20
volume = 30*20*10
CUBE = [x,y,z]
STEEL = {
    "yield": 350,
    "elastic": 210000,
    "shear":81000,
    "poisson":0.3
}
ALUMINIUM = {
    "yield": 270,
    "elastic": 70000,
    "shear":27000,
    "poisson":0.3
}
# data points: 
points = [[0,0,20],[30,0,20],[30,10,20],[0,10,20],[0,0,0],[30,10,0],[0,10,0],[0,5,0],[0,5,20],[30,5,20],[30,5,0]]
# need to respect order (same index in points and displacements)
after_loading = [[0.001,0.002,20],[30.001,0,20.004],[29.997,10.003,19.996],[0.004,10.009,19.995],[0,0,0],[30.003,0.001,0.0026],[29.996,10.003,0],[0.001,9.996,0.0021],[0,4.995,0],[0.004,5.009,19.996],[29.997,5.003,19.996],[29.996,5.003,0]] #[u,v,w]
def strain_transform(point,afterLoad,MATERIAL) :
    displacement=[]
    for x in range(len(point)):
        displacement.append(afterLoad[x]-point[x])
    c = [0,0,0]
    eps = [0,0,0] 
    for i in range(3):
        c[i] = (displacement[i] - point[i]) / volume
        eps[i] = displacement[i]/CUBE[i]
    e = sum(eps)
    gxy = c[0]*x*z + c[1]*y*z
    gxz = c[0]*x*y + c[2]*y*z 
    gyz = c[2]*x*z + c[1]*x*y
    d =[0,0,0]
    for i in range(len(eps)):
        d[i] = 2*MATERIAL["shear"]*eps[i] + get_lambda(MATERIAL)*e
    txy = MATERIAL["shear"]*2*gxy
    txz = MATERIAL["shear"]*2*gxz
    tyz = MATERIAL["shear"]*2*gyz
    strain_matrix = [ 
        [eps[0],gxy,gxz],
        [gxy,eps[1],gyz],
        [gxz,gyz,eps[2]]
    ]
    stress_matrix = [
        [d[0],txy,txz],
        [txy,d[1],tyz],
        [txz,tyz,d[2]]
        ]
    oct_stress = (d[0] - d[1])**2 + (d[1]-d[2])**2 + 6*(txy**2 + txz**2 +tyz**2)
    oct_stress = sqrt(oct_stress) / 3
    uod = (3 * oct_stress**2) / (4*MATERIAL["shear"])
    u0 = sum([d[0]*eps[0], d[1] * eps[1], d[2] * eps[2],txy*gxy,txz*gxz,tyz*gyz]) /2
    v,d= la.eig(stress_matrix)
    return {
       "strain" : strain_matrix,
        "stess": stress_matrix,
        "v":v,
        "d":d,
        "octahedral_stress":oct_stress,
        "distortional_strain":uod,
        "energy_density":u0
    }

def get_lambda(MATERIAL):
    return (MATERIAL["poisson"]*MATERIAL["elastic"]/((1+MATERIAL["poisson"])*(1-2*MATERIAL["poisson"])))

def save(filename, results):
    print(filename)
    with open(filename, 'w',newline='') as file:
        file.write('-'*30 + '\n')
        file.write('{:<17}'.format("name") + " | " + "value\n")
        for (key,value) in results.items():
            file.write('\n{:<17}'.format(key) + " |\n")
            file.write(str(value))
            file.write('\n')
            file.write("-"*30)

if __name__ == "__main__":
    for x in range (len(points)):
        if x in (2,3,6):
            material = STEEL
        else :
            material = ALUMINIUM
        result = strain_transform(points[x],after_loading[x], material)
        save(f"point_{x}.txt", result)
    for x in range (len(points)):
        if x in (8,9,10,11,12,13):
            material = STEEL
        result = strain_transform(points[x],after_loading[x], material)
        save(f"point[mid]_{x}.txt", result)