"""
    microphase(::RootsOptimizer, m::MicrophaseModel)

Compute the microphase separation boundary between two phases. The phase boundary between these two phases is not necessary the true phase boundary because other phases may have lower free energy than current two phases. Therefore, it relies on the user to make sure whether these two phases are actually two lowest phases in this parameter region.
"""
function microphase(ro::RootsOptimizer, m::MicrophaseModel)
    F1 = x -> free_energy(x, 1, m)[1]
    F2 = x -> free_energy(x, 2, m)[1]
    tracker = Roots.Tracks()
    xstar = Roots.find_zero(ϕ -> F1(ϕ) - F2(ϕ), (m.a, m.b), ro.algo;
                            tracks=tracker, ro.kwargs...)
    return xstar, tracker
end

"""
    microphase(::RootsOptimizer, model::SurrogateMicrophaseModel)

Compute the microphase separation boundary between two phases whose free energies are surrogated by the input `SurrogateMicrophaseModel`.

Inside `SurrogateMicrophaseModel` model, use `CHSplineSurrogate` for `ϕParam` and `SplineSurrogate` for other type of `PolymerParameter`s.
"""
function microphase(ro::RootsOptimizer, model::SurrogateMicrophaseModel)
    s1 = model.surrogate1
    s2 = model.surrogate2
    F1 = ϕ -> free_energy(ϕ, s1)[1]
    F2 = ϕ -> free_energy(ϕ, s2)[1]
    # Find the common range for both phases which make sure the surrogates work.
    ϕa = max(s1.x[1], s2.x[1])
    ϕb = min(s1.x[end], s2.x[end])
    tracker = Roots.Tracks()
    xstar = Roots.find_zero(ϕ -> F1(ϕ) - F2(ϕ), (ϕa, ϕb), ro.algo;
                            tracks=tracker, ro.kwargs...)
    return xstar, tracker
end

"""
    microphase(::AbstractSurrogateSolver, m::MicrophaseModel)

Compute the microphase separation boundary between two phases using AbstractSurrogateSolver. The phase boundary between these two phases is not necessary the true phase boundary because other phases may have lower free energy than current two phases. Therefore, it relies on the user to make sure whether these two phases are actually two lowest phases in this parameter region.

IMPORTANT NOTE: due to lack of particular accurate surrogate model for microphase process (unlike macrophase case), the converence rate is typically slow than RootsOptimizer.
"""
function microphase(os::AbstractSurrogateSolver, m::MicrophaseModel)
    F1 = x -> free_energy(x, 1, m)[1]
    F2 = x -> free_energy(x, 2, m)[1]
    xstar, tracker = microphase(os, m.ctrlparam, F1, F2)
    return xstar, tracker
end

function microphase(os::OPSOptimizer, param, f1, f2)
    vecx1 = deepcopy(os.x01)
    vecx2 = deepcopy(os.x02)
    if isnothing(os.f01)
        vecfg1 = f1.(vecx1)
        vecfg2 = f2.(vecx2)
        vecf1 = first.(vecfg1)
        vecf2 = first.(vecfg2)
    else
        vecf1 = deepcopy(os.f01)
        vecf2 = deepcopy(os.f02)
    end
    vecxans = copy(vecx1)
    vecf1ans = copy(vecf1)
    vecf2ans = copy(vecf2)

    for i in 1:os.max_itr
        surrogate1 = SplineSurrogate(vecx1, vecf1; k=2)
        surrogate2 = SplineSurrogate(vecx2, vecf2; k=2)

        smpmodel = SurrogateMicrophaseModel(surrogate1, surrogate2, param)

        x, _ = microphase(os.optimizer, smpmodel)

        push!(vecxans, x)

        fn1 = f1(x)
        push!(vecf1ans, fn1)
        push!(vecx1, x)
        push!(vecf1, fn1)
        p1 = sortperm(vecx1)
        vecx1, vecf1 = vecx1[p1], vecf1[p1]

        fn2 = f2(x)
        push!(vecf2ans, fn2)
        push!(vecx2, x)
        push!(vecf2, fn2)
        p2 = sortperm(vecx2)
        vecx2, vecf2 = vecx2[p2], vecf2[p2]

        if i > 1 && abs(x - vecxans[end-1]) < os.xatol
            break
        end
    end

    niter = length(vecxans) - length(os.x01)
    return vecxans[end], (2*niter, vecxans, vecf1ans, vecf2ans)
end