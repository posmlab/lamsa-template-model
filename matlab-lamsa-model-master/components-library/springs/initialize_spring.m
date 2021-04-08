%% function to initialize a new spring object

function spring = initialize_spring(app)

    type = app.spring.SelectedTab;
    
    if (type == app.linear_spring)
        spring = LinearSpring(app.linear_spring_k.Value,app.linear_spring_mass.Value,app.linear_spring_Fmax.Value);
    elseif (type == app.exponential_spring)
        spring = ExponentialSpring(app.exp_spring_k.Value,app.exp_spring_char_len.Value,app.exp_spring_mass.Value,app.exp_spring_Fmax.Value);
    elseif (type == app.linear_elastic_extensional_spring)
        spring = LinearElasticExtensionalSpring(app.lee_spring_E.Value,app.lee_spring_A.Value,app.lee_spring_L.Value,app.lee_spring_rho.Value,app.lee_spring_sigma_f.Value);
    end

end