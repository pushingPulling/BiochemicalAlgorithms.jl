function train(; cuda = true, η₀ = 1.0f-3, λ = 1.0f-4, epochs = 500)
    if cuda && CUDA.has_cuda()
        device = gpu
        CUDA.allowscalar(false)
        @info "Training on GPU"
    else
        device = cpu
        @info "Training on CPU"
    end

    model = FourierNeuralOperator(ch = (2, 64, 64, 64, 64, 64, 128, 1), modes = (16,),
                                  σ = gelu)
    data = get_2d_dataloader()
    optimiser = Flux.Optimiser(WeightDecay(λ), Flux.Adam(η₀))
    loss_func = l₂loss

    learner = Learner(model, data, optimiser, loss_func,
                      ToDevice(device, device))

    fit!(learner, epochs)
    model = learner.model |> cpu
    @save "model/model_burger.bson" model

    return learner
end

function get_2d_dataloader()
    train_set = Dataloader()
end