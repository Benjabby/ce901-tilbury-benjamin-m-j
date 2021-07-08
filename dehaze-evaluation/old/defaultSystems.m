function systems = defaultSystems()

    systems = struct();
    
    %systems.Ancuti = @dehazeAncuti;
    systems.Berman = @dehazeBerman;
    systems.Cai = @dehazeCai;
    systems.He = @dehazeHe;
    systems.Tarel = @dehazeTarel;
    systems.Tilbury = @dehazeTilbury;
    systems.Tsai = @dehazeTsai;
    systems.Zhu = @dehazeZhu;
    
    systems = prepareSystems(systems);
end