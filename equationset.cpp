#include <future>
#include <thread>
#include "equationset.h"
#include "physics_cubic.h"
#include "physics_hex.h"
#include <chrono>
#include <iostream>

EquationSet::EquationSet(bool cubic)
{
    this->cubic = cubic;
}


template<class F>
struct Forward {
    F f;
    template<class... Ts>
    auto operator()(Ts&&... as) {
        //double ratio = std::chrono::high_resolution_clock::period::num/double(std::chrono::high_resolution_clock::period::den);
        //auto start = std::chrono::high_resolution_clock::now();
        auto res = f(as...);
        //auto end = std::chrono::high_resolution_clock::now();
        //std::cout<<(end-start).count()*ratio*1000<<std::endl;
        return res;
    }
};

template<class F>
auto timed(F f) {
    return Forward<F> {f};
}

RealVector EquationSet::evaluate(RealVector point, double L){
    if(cubic){
        double zp_der = Zp_der_cubic(point, L);
/*
        std::packaged_task<double(RealVector, double, double)> a_task(timed(A_der_cubic)),
                ms_task(timed(Ms_der_cubic)), mp_task(timed(Mp_der_cubic)), zs_task(timed(Zs_der_cubic));
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
                                              zs_der.get(), zp_der, 0});
*/

        return RealVector(std::vector<double>{A_der_cubic(point, L, zp_der),
                                              Ms_der_cubic(point, L, zp_der),
                                              Mp_der_cubic(point, L, zp_der),
                                              Zs_der_cubic(point, L, zp_der),
                                              zp_der,
                                              0.});

    } else {
        double zp_der = Zp_der_hex(point, L);
        return RealVector(std::vector<double>{A_der_hex(point, L, zp_der),
                                              Ms_der_hex(point, L, zp_der),
                                              Mp_der_hex(point, L, zp_der),
                                              Zs_der_hex(point, L, zp_der),
                                              zp_der,
                                              0.});
    }
}
