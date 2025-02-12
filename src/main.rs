use std::io;
use std::f64::consts::PI;
#[derive(Debug)]
struct Vector3D {
    x:f64,
    y:f64,
    z:f64,
}
struct PointCurvilinear{
    fi:f64,
    la:f64,
    h:f64,
}

impl Vector3D {
    fn add(&self, other:&Vector3D)-> Vector3D{
        Vector3D {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
    fn subtract (&self, other:&Vector3D)-> Vector3D{
        Vector3D {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
    fn module(&self) -> f64 {
        ((self.x * self.x) + (self.y * self.y) + (self.z * self.z)).sqrt()
    }
}




fn read_vector(prompt: &str) -> Vector3D {
    println!("{}", prompt);
    let mut input = String::new();

    // Leer la componente x
    println!("Ingrese la componente x:");
    io::stdin().read_line(&mut input).expect("Error al leer la entrada");
    let x: f64 = input.trim().parse().expect("La componente x debe ser un número");
    input.clear(); // Limpiar el buffer de entrada

    // Leer la componente y
    println!("Ingrese la componente y:");
    io::stdin().read_line(&mut input).expect("Error al leer la entrada");
    let y: f64 = input.trim().parse().expect("La componente y debe ser un número");
    input.clear(); // Limpiar el buffer de entrada

    // Leer la componente z
    println!("Ingrese la componente z:");
    io::stdin().read_line(&mut input).expect("Error al leer la entrada");
    let z: f64 = input.trim().parse().expect("La componente z debe ser un número");

    Vector3D { x, y, z }
}
fn read_coordinate(prompt: &str) -> PointCurvilinear {
    println!("{}", prompt);
    let mut input = String::new();

    // Leer la componente x
    println!("Ingrese la componente latitud en grados decimales:");
    io::stdin().read_line(&mut input).expect("Error al leer la entrada");
    let fi: f64 = input.trim().parse().expect("La componente fi debe ser un número");
    input.clear(); // Limpiar el buffer de entrada

    // Leer la componente y
    println!("Ingrese la componente longitud en grados decimales:");
    io::stdin().read_line(&mut input).expect("Error al leer la entrada");
    let la: f64 = input.trim().parse().expect("La componente lambda debe ser un número");
    input.clear(); // Limpiar el buffer de entrada

    // Leer la componente z
    println!("Ingrese la componente altura elipsoidal en metros:");
    io::stdin().read_line(&mut input).expect("Error al leer la entrada");
    let h: f64 = input.trim().parse().expect("La componente h debe ser un número");

    PointCurvilinear { fi, la, h }
}
fn computeN (a:f64, b:f64, punto:&PointCurvilinear)->f64{
    let fi_radianes = punto.fi * (PI / 180.0); // Convertir grados a radianes
    let numerador = a * a;
    println!("numerador es: {}", numerador);
    
    let denominador = (((a * a) * fi_radianes.cos().powi(2)) + ((b * b) * fi_radianes.sin().powi(2))).sqrt();
    println!("numerador es: {}", denominador);
    numerador / denominador
}
fn computeCartesian(punto:&PointCurvilinear,primerVertical:f64, a:f64, b:f64)->Vector3D{
    let fi_radianes = punto.fi * (PI / 180.0);
    let lambda_radianes = punto.la * (PI / 180.0);
    let x = (primerVertical+punto.h)* fi_radianes.cos()*lambda_radianes.cos();
    let y = (primerVertical+punto.h)* fi_radianes.cos()*lambda_radianes.sin();
    let z = ((primerVertical * b.powi(2) / a.powi(2))+punto.h)* fi_radianes.sin();
    Vector3D {x,y,z}
}

fn main() {
    //WGS84 parameters
    let a = 6378137.000;
    let b = 6356752.3142;
    // Pedir al usuario los dos vectores
    //let v1 = read_vector("Ingrese el primer vector:");
    //let v2 = read_vector("Ingrese el segundo vector:");
    //Pedir al usuario dos coordenadas 
    let p1  =read_coordinate("Ingrese las coordenadas del primer punto");
    let p2  =read_coordinate("Ingrese las coordenadas del segundo punto");
    //Calcular los vectores Normal de cada punto 
    let n1 = computeN(a,b,&p1);
    let n2 = computeN(a,b,&p2);
    //Obtener las coordenadas cartesianas de cada punto
    let v1 = computeCartesian(&p1,n1,a,b);
    let v2 = computeCartesian(&p2,n2,a,b);
    let resultado_resta = v1.subtract(&v2);
    let modulo = resultado_resta.module();
    println!("La primer vertical n1 es: {}", n1);
    println!("La segunda vertical n1 es: {}", n2);
    println!("Las coordenadas cartesianas de la cdmx son: {:?}", v1);
    println!("Las coordenadas cartesianas de toluca son: {:?}", v2);
    println!("El vector formado de ambos puntos es: {:?}", resultado_resta);
    println!("El módulo del vector resultante es: {}", modulo);

}