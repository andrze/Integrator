#include <future>
#include <thread>
#include "equationset.h"
#include "physics.h"


EquationSet::EquationSet()
{

}

RealVector EquationSet::evaluate(RealVector point, double L){

    double zp_der = Zp_der(point, L);
/*
    std::packaged_task<double(RealVector, double, double)> a_task(A_der),
            ms_task(Ms_der), mp_task(Mp_der), zs_task(Zs_der);
    std::future<double> a_der=a_task.get_future(),
            ms_der=ms_task.get_future(), mp_der=mp_task.get_future(),
            zs_der=zs_task.get_future();
    std::thread t_a(std::move(a_task), point, L, zp_der),
            t_ms(std::move(ms_task), point, L, zp_der),
            t_mp(std::move(mp_task), point, L, zp_der),
            t_zs(std::move(zs_task), point, L, zp_der);
    t_a.join();
    t_ms.join();
    t_mp.join();
    t_zs.join();
    return RealVector(std::vector<double>{a_der.get(),
                                          ms_der.get(), mp_der.get(),
                                          zs_der.get(), zp_der});

    double a_der = A_der(point, L, zp_der);
    double ms_der = Ms_der(point, L, zp_der);
    double mp_der = Mp_der(point, L, zp_der);
    double zs_der = Zs_der(point, L, zp_der);
*/
    return RealVector(std::vector<double>{A_der(point, L, zp_der),
                                          Ms_der(point, L, zp_der),
                                          Mp_der(point, L, zp_der),
                                          Zs_der(point, L, zp_der),
                                          zp_der,
                                          0.});
}
