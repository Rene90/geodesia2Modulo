use std::io;
use std::f64::consts::PI;
#[derive(Debug)]
struct Vector3D {
    x:f64,
    y:f64,
    z:f64,
}
struct Datadist {
    azimuth12:f64,
    dist12:f64,
}
struct PointCurvilinear{
    fi:f64,
    la:f64,
    h:f64,
}
struct PointCurvilinear1{
    h:f64,
    phi:f64,
    n:f64,
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


fn read_azdist(prompt: &str) -> Datadist {
    println!("{}", prompt);
    let mut input = String::new();
    println!("Ingrese la distancia al punto 2:");
    io::stdin().read_line(&mut input).expect("Error al leer la entrada");
    let dist12: f64 = input.trim().parse().expect("La distancia debe ser un número");
    input.clear(); // Limpiar el buffer de entrada
    println!("Ingrese el azimuth al punto 2:");
    io::stdin().read_line(&mut input).expect("Error al leer la entrada");
    let azimuth12: f64 = input.trim().parse().expect("el azimuth debe ser un numero");
    input.clear(); // Limpiar el buffer de entrada
    Datadist{ azimuth12,dist12 }
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
fn compute_n(a: f64, b: f64, p_geo: &PointCurvilinear1) -> f64 {
    let e2 = 1.0 - (b.powi(2) / a.powi(2));
    a / (1.0 - e2 * p_geo.phi.sin().powi(2)).sqrt()
}
fn compute_nsencillo(a: f64, b: f64, latitud: f64) -> f64 {
    let e2 = 1.0 - (b.powi(2) / a.powi(2));
    a / (1.0 - e2 * latitud.sin().powi(2)).sqrt()
}
fn compute_msencillo(a: f64, b: f64, latitud:f64) -> f64 {
    let e2 = 1.0 - (b.powi(2) / a.powi(2));
    (a*(1.0 - e2)) / (1.0 - e2 * latitud.sin().powi(2)).powf(1.5)
}
fn compute_m(a: f64, b: f64, p_geo: &PointCurvilinear) -> f64 {
    let e2 = 1.0 - (b.powi(2) / a.powi(2));
    (a*(1.0 - e2)) / (1.0 - e2 * p_geo.fi.sin().powi(2)).powf(1.5)
}
fn computeN (a:f64, b:f64, punto:&PointCurvilinear)->f64{
    let fi_radianes = punto.fi * (PI / 180.0); // Convertir grados a radianes
    let numerador = a * a;
    println!("numerador es: {}", numerador);
    
    let denominador = (((a * a) * fi_radianes.cos().powi(2)) + ((b * b) * fi_radianes.sin().powi(2))).sqrt();
    println!("denominador es: {}", denominador);
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
fn computeLongitude(punto: &Vector3D)-> f64{
    let y = punto.y;
    let x = punto.x;
    let arctan2_yx = y.atan2(x);
    arctan2_yx 
}

fn compute_iterative(punto: &Vector3D, a: f64, b: f64) -> (f64, f64, f64) {
    // Valores iniciales
    let mut n = a;
    let z = punto.z;
    let x = punto.x;
    let y = punto.y;
    let e2 = 1.0 - (b.powi(2) / a.powi(2));
    let p = (x.powi(2) + y.powi(2)).sqrt();
    let magnitud = punto.module();
    let mut h = magnitud - (a * b).sqrt();
    let phi1 = z / p;
    let mut phi2 = 1.0 - (e2 * n / (n + h));
    let mut phi3 = phi1 / phi2;
    let mut phi = phi3.atan();
    let mut contador =0;

    // Umbral de convergencia
    let threshold = 1e-6;

    // Bucle iterativo
    loop {
        contador = contador + 1;
        // Guardar los valores anteriores
        let h_prev = h;
        let phi_prev = phi;

        // Calcular nuevos valores
        n = compute_n(a, b, &PointCurvilinear1 { h, phi, n });
        h = compute_h_iterative(punto, &PointCurvilinear1 { h, phi, n }, phi, a, b);
        phi = compute_latitude_iterative(punto, &PointCurvilinear1 { h, phi, n }, h, a, b);
        println!("iteracion: {} ", contador);
        println!("phi: {} grados", phi * 180.0 / PI);
        println!("h: {} ", h);
        println!("N: {} ", n);
        // Verificar convergencia
        if (h - h_prev).abs() < threshold && (phi - phi_prev).abs() < threshold {
            break;
        }
    }

    (h, phi, n)
}

// Función para calcular h iterativamente
fn compute_h_iterative(punto: &Vector3D, p_geo: &PointCurvilinear1, phi: f64, a: f64, b: f64) -> f64 {
    let z = punto.z;
    let x = punto.x;
    let y = punto.y;
    let fphi = phi; // phi ya está en radianes
    let p = (x.powi(2) + y.powi(2)).sqrt();
    let n = compute_n(a, b, p_geo);
    let h1 = p / fphi.cos();
    h1 - n
}

// Función para calcular la latitud iterativamente
fn compute_latitude_iterative(punto: &Vector3D, p_geo: &PointCurvilinear1, h: f64, a: f64, b: f64) -> f64 {
    let z = punto.z;
    let x = punto.x;
    let y = punto.y;
    let e2 = 1.0 - (b.powi(2) / a.powi(2));
    let p = (x.powi(2) + y.powi(2)).sqrt();
    let n = compute_n(a, b, p_geo);
    let z_p = z / p;
    let temp_e = 1.0 - (e2 * n / (n + h));
    let temp_inv = 1.0 / temp_e;
    let ze = z_p * temp_inv;
    ze.atan()
}
//Funcion para calcular el arco de meridiano entre 2 latitudes diferentes usando aproximacion por expansion de Binomio de Newton y solucionando de manera analitica
fn arcoMeridiano(p_geo1: &PointCurvilinear, p_geo2: &PointCurvilinear, a: f64, b:f64)->f64{

    let phi1_rad = p_geo1.fi* (PI / 180.0);
    let phi2_rad = p_geo2.fi* (PI / 180.0);
    let delta_phi = phi2_rad - phi1_rad;
    let e2 = 1.0 - (b.powi(2) / a.powi(2));

    // Términos de la expansión en serie
    let term1 = delta_phi;
    let term2 = (3.0 / 4.0) * e2 * delta_phi;
    let term3 = -(3.0 / 8.0) * e2 * (f64::sin(2.0 * phi2_rad) - f64::sin(2.0 * phi1_rad));
    let term4 = (15.0 / 64.0) * e2.powi(2) * delta_phi;
    let term5 = -(15.0 / 32.0) * e2.powi(2) * (f64::sin(2.0 * phi2_rad) - f64::sin(2.0 * phi1_rad));
    let term6 = (15.0 / 256.0) * e2.powi(2) * (f64::sin(4.0 * phi2_rad) - f64::sin(4.0 * phi1_rad));

    // Sumar todos los términos
    let s = a * (1.0 - e2) * (term1 + term2 + term3 + term4 + term5 + term6);
    s
}
fn arcoMeridianoPromedio(p_geo1: &PointCurvilinear, p_geo2: &PointCurvilinear, a: f64, b:f64)->f64{

    let phi1_rad = p_geo1.fi* (PI / 180.0);
    let phi2_rad = p_geo2.fi* (PI / 180.0);
    let phimedia = (phi1_rad+phi2_rad)/2.0;
    let delta_phi = phi2_rad - phi1_rad;
    let e2 = 1.0 - (b.powi(2) / a.powi(2));
    let m = compute_msencillo(a,b,phimedia);
    let s= delta_phi*m;
    s


   
}

//Funcion para calcular el arco de paralelo ente 2 longitudes diferentes radio de circulo utilizando la primer vertical o N
fn arcoParalelo(p_geo1: &PointCurvilinear, p_geo2: &PointCurvilinear,a: f64, b:f64)->f64{
    let phi1_rad = p_geo1.fi* (PI / 180.0);
    let phi2_rad = p_geo2.fi* (PI / 180.0);
    let la1_rad = p_geo1.la * (PI / 180.0);
    let la2_rad = p_geo2.la * (PI / 180.0);
    let menor = phi1_rad.min(phi2_rad);
    //Utilizo la latiud menor para sacar el radio de la primer normal porque es el cateto adyacente del triangulo rectangulo formado entre ambos puntos
    let nm = compute_nsencillo(a,b,menor);
    println!("nm {}", nm);
    let radio = nm*menor.cos();
    //diferencia de longitudes
    let deltaLambda=  la2_rad-la1_rad;
    let s = radio * deltaLambda;
    s
}

//Funcio para calcular distancia entre dos puntos usando el arco de paralelo y el arco de meridiano como los catetos de un  triangulo rectangulo
fn calcularDistanciaArcos(arcoM: f64, arcoP: f64)-> f64{
    let s = arcoM.powi(2)+ arcoP.powi(2);
    let s2 =s.sqrt();
    s2
}
// Funcion para calcular distancia y azimuth con Puissant a traves de iteraciones
fn computePuissantInverso(p_geo1: &PointCurvilinear,p_geo2: &PointCurvilinear, a: f64, b: f64)-> (f64, f64) {
    let n1 = computeN(a, b, p_geo1);
    let m1 = compute_m(a, b, p_geo1);
    let n2 = computeN(a, b, p_geo2);
    let phi1=p_geo1.fi* (PI / 180.0);
    let phi2=p_geo2.fi* (PI / 180.0);
    let lam1=p_geo1.la* (PI / 180.0);
    let lam2=p_geo2.la* (PI / 180.0);
    let delta_lambda = (lam1 - lam2).abs();
    let delta_phi = (phi1 - phi2).abs();
    let e2 = 1.0 - (b.powi(2) / a.powi(2));
    // primera iteracion, se usa primer termino para aproximar el azimuth directo
    let primeralpha1216= delta_phi * m1 / (1.0 - (3.0 * e2 * phi1.sin() * phi1.cos()) / (2.0 * (1.0 - e2 * phi1.sin().powi(2))));    
    let complemento_phi = (PI/2.0)-phi2;
    let primeralpha1217 =delta_lambda*n2*complemento_phi.sin();
    let mut cocientePrimero = (primeralpha1217 / primeralpha1216).atan();
    // primera distancia
    let mut primers1217 = delta_lambda*n2*complemento_phi.sin()/cocientePrimero;
    println!("Azimuth primera aprox : {}", cocientePrimero* (180.0 / PI));
    println!("Distancia primera aprox : {}", primers1217);  
    // iteraciones
    let thresholdAzimuth = 1e-6;
    let thresholdDistancia = 1e-3;
    let mut iter_count = 0;
    let max_iter = 15;
    loop {
        iter_count += 1;
        if iter_count > max_iter {
            println!("Advertencia: No converge después de {} iteraciones", max_iter);
            break;
        }
        let term2 = primers1217.powi(2) * phi1.tan() * primeralpha1217.sin().powi(2) / (2.0 * m1 * n1);
        let term3 = primers1217.powi(3) * primeralpha1217.cos() * primeralpha1217.sin().powi(2) * (1.0 + 3.0 * phi1.tan().powi(2)) / (6.0 * m1 * n1.powi(2));
        let alpha1216 = ((delta_phi / (1.0 - (3.0 * e2 * phi1.sin() * phi1.cos()) / (2.0 * (1.0 - e2 * phi1.sin().powi(2)))))+ term2 + term3) * m1;
        let term2l = (primers1217.powi(3)/(6.0*n2.powi(3))) * (primeralpha1217.sin()/complemento_phi.sin());
        let term3l = (primers1217.powi(3)/(6.0*n2.powi(3))) * (primeralpha1217.sin().powi(3)/complemento_phi.sin().powi(3));
        let alpha1217 = (delta_lambda + term2l - term3l)*n2*complemento_phi.sin();
        let cocientesegundo = (alpha1217 / alpha1216).atan();
        let segundos1217  = (delta_lambda + term2l - term3l)*n2 /  (cocientesegundo.sin()/complemento_phi.sin());
        println!("Iteración: Azimuth: {}, Distancia: {}", cocientesegundo * (180.0 / PI), segundos1217);

        if (cocientesegundo - cocientePrimero).abs() < thresholdAzimuth && (segundos1217 - primers1217).abs()  < thresholdDistancia{
            cocientePrimero = cocientesegundo;
            primers1217 =segundos1217;
            break;
        }
        cocientePrimero = cocientesegundo;
        primers1217 =segundos1217;   

    }

    (cocientePrimero,primers1217)


}

// Funcion para calcular el punto siguiente con la formula de Puissant hasta 150 kilometros
fn compute_p2_Puissant(p_geo: &PointCurvilinear, a: f64, b: f64, s: f64, azimuth: f64)-> PointCurvilinear {
    let n = computeN(a, b, p_geo);
    let m = compute_m(a, b, p_geo);
    let e2 = 1.0 - (b.powi(2) / a.powi(2));
    let phi1=p_geo.fi* (PI / 180.0);
    let lam1=p_geo.la* (PI / 180.0);
    let s12=s;
    let alpha12 = azimuth* (PI / 180.0);
    //Primero calculamos delta phi para obtener phi2
    let term1 = s12 * alpha12.cos() / m;
    let term2 = s12.powi(2) * phi1.tan() * alpha12.sin().powi(2) / (2.0 * m * n);
    let term3 = s12.powi(3) * alpha12.cos() * alpha12.sin().powi(2) * (1.0 + 3.0 * phi1.tan().powi(2)) / (6.0 * m * n.powi(2));
    let rho_prime_prime = term1 - term2 - term3;
    //let s1216 = ((delta_phi / (1.0 - (3.0 * e2 * phi1.sin() * phi1.cos()) / (2.0 * (1.0 - e2 * phi1.sin().powi(2))))) + term2 + term3)* m/ alpha12.cos();
    
    let delta_phi = rho_prime_prime * (1.0 - (3.0 * e2 * phi1.sin() * phi1.cos()) / (2.0 * (1.0 - e2 * phi1.sin().powi(2))));
    //let alpha1216 = ((delta_phi / (1.0 - (3.0 * e2 * phi1.sin() * phi1.cos()) / (2.0 * (1.0 - e2 * phi1.sin().powi(2)))))+ term2 + term3) * m;
    //let primeralpha1216= delta_phi * m / (1.0 - (3.0 * e2 * phi1.sin() * phi1.cos()) / (2.0 * (1.0 - e2 * phi1.sin().powi(2))));
    let phi2 = delta_phi* (180.0 / PI) + p_geo.fi;
    println!("La delta Phi es de : {}", delta_phi* (180.0 / PI));
    //calculamos delta lambda para obtener lambda2
    let latitud_media = (p_geo.fi+phi2)/2.0;
    let nm = compute_nsencillo(a,b,phi2);
    let complemento_phi = (PI/2.0)-phi1;
    let term1l = (s12/nm) * (alpha12.sin()/complemento_phi.sin());
    let term2l = (s12.powi(3)/(6.0*nm.powi(3))) * (alpha12.sin()/complemento_phi.sin());
    let term3l = (s12.powi(3)/(6.0*nm.powi(3))) * (alpha12.sin().powi(3)/complemento_phi.sin().powi(3));
    
    //let s1217  = (delta_lambda + term2l - term3l)*nm /  (alpha12.sin()/complemento_phi.sin());
    //let alpha1217 = (delta_lambda + term2l - term3l)*nm /  (1/complemento_phi.sin());
    //let primeralpha1217 =delta_lambda*nm*complemento_phi.sin();
    //let primers1217 = delta_lambda*nm*complemento_phi.sin()/primeralpha1217;


    let delta_lambda = term1l-term2l+term3l;
    println!("La delta Lambda es de : {}", delta_lambda* (180.0 / PI));
    let lambda2 = delta_lambda* (180.0 / PI) + p_geo.la;
    let height =0.0;
    PointCurvilinear{fi:phi2,la:lambda2,h:p_geo.h}
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
    // let p2  =read_coordinate("Ingrese las coordenadas del segundo punto");


    //let (resultadoAzimuth,resultadodistancia)=computePuissantInverso(&p1,&p2, a, b);

    //problema directo punto , distancia y azimuth
    let azdist = read_azdist("Ingrese la distancia y azimuth al siguiente punto");
    let azimuth=azdist.azimuth12;
    let s= azdist.dist12;
    //let p2  =read_coordinate("Ingrese las coordenadas del segundo punto");
    //Calcular los vectores Normal de cada punto 
    let n1 = computeN(a,b,&p1);
    let p2 =compute_p2_Puissant(&p1, a, b, s, azimuth);
    println!("La latitud del punto 2  es: {}", p2.fi);
    println!("La longitud del punto 2  es: {}", p2.la);
    //Calcular la distancia entre ambos puntos para comprobar aproximacion, usando los arcos de meridiano y de paralelo
    //let arcoMer =arcoMeridiano(&p1, &p2, a, b);
    //let arcoMer =arcoMeridianoPromedio(&p1, &p2, a, b);
    //let arcoPar =arcoParalelo(&p1, &p2, a, b);
    //let distanciaPuntos = calcularDistanciaArcos(arcoMer,arcoPar);
    //println!("El arco de meridiano entre ambos puntos es de : {}", arcoMer);
    //println!("El arco de paralelo entre ambos puntos es de {}", arcoPar);
    //println!("La distancia entre ambos puntos usando los arcos de meridiano y paralelo es de  {}", distanciaPuntos);
    //println!("La diferencia entre la distancia original y la calculada es de {}", s-distanciaPuntos);
    //let n2 = computeN(a,b,&p2);
    //let v1 = computeCartesian(&p1,n1,a,b);
    //let v2 = computeCartesian(&p2,n2,a,b);
    //let resultado_resta = v1.subtract(&v2);
    //let modulo = resultado_resta.module();
    //println!("La distancia entre ambos puntos usando la distancia euclideana es de  {}", modulo);
    //println!("La diferencia entre la distancia original y la geometrica es de {}", s-modulo);



    //let n2 = computeN(a,b,&p2);
    //Obtener las coordenadas cartesianas de cada punto
    //let v1 = computeCartesian(&p1,n1,a,b);
    //let v2 = computeCartesian(&p2,n2,a,b);
    //let resultado_resta = v1.subtract(&v2);
    //let modulo = resultado_resta.module();
    //println!("La primer vertical n1 es: {}", n1);
    //println!("La segunda vertical n1 es: {}", n2);
    //println!("Las coordenadas cartesianas del primer punto son: {:?}", v1);
    //println!("Las coordenadas cartesianas del segundo punto son: {:?}", v2);
    //println!("El vector formado de ambos puntos es: {:?}", resultado_resta);
    //println!("El módulo del vector resultante es: {}", modulo);
    //println!("Calculando las coordenadas geodesicas de vuelta");
    //let longitud1 = computeLongitude(&v1);
    //println!("longitud1: {}", longitud1);
    
    //let longitud2 =computeLongitude(&v2);
    //println!("longitud2: {}", longitud2);
    //let (h11, phi11, n11) = compute_iterative(&v1, a, b);
    //println!("h: {}", h11);
    //println!("phi: {} radianes", phi11);
    //println!("phi: {} grados", phi11 * 180.0 / PI);
    //println!("n: {}", n11);
    //let (h22, phi22, n22) = compute_iterative(&v2, a, b);
    //println!("h: {}", h22);
    //println!("phi: {} radianes", phi22);
    //println!("phi: {} grados", phi22 * 180.0 / PI);
    //println!("n: {}", n22);



}